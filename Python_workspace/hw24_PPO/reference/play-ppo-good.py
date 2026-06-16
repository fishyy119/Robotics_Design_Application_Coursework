# deploy_ppo_render.py
# -*- coding: utf-8 -*-
"""
根据 ppo-good.py 训练得到的 ppo_model.pt 进行部署/测试，并使用 render 显示环境。

默认适配：
    CartPole-v1 + ppo_model.pt

运行示例：
    python deploy_ppo_render.py --env CartPole-v1 --model ppo_model.pt --episodes 5 --device cpu

如果训练时使用了像素输入：
    python deploy_ppo_render.py --env Pong-v0 --model ppo_model.pt --from_pixel --episodes 3
"""

import argparse
import time
import numpy as np
import gymnasium as gym
import torch
import torch.nn as nn
from gymnasium.spaces import Box, Discrete
from torch.distributions.normal import Normal
from torch.distributions.categorical import Categorical


# =========================
# 下面这些类名和结构需要与训练脚本保持一致
# 因为训练脚本使用 torch.save(ac, 'ppo_model.pt') 保存的是完整对象
# =========================

def cnn(in_channels, activation=nn.ReLU):
    return nn.Sequential(
        nn.Conv2d(in_channels=in_channels, out_channels=16, kernel_size=8, stride=4),
        activation(),
        nn.Conv2d(in_channels=16, out_channels=32, kernel_size=4, stride=2),
        activation(),
        nn.Flatten()
    )


def mlp(sizes, activation, output_activation=nn.Identity):
    layers = []
    for j in range(len(sizes) - 1):
        act = activation if j < len(sizes) - 2 else output_activation
        layers += [nn.Linear(sizes[j], sizes[j + 1]), act()]
    return nn.Sequential(*layers)


class Actor(nn.Module):
    def _distribution(self, obs):
        raise NotImplementedError

    def _log_prob_from_distribution(self, pi, act):
        raise NotImplementedError

    def forward(self, obs, act=None):
        pi = self._distribution(obs)
        logp_a = None
        if act is not None:
            logp_a = self._log_prob_from_distribution(pi, act)
        return pi, logp_a


class CategoricalActor(Actor):
    def __init__(self, obs_dim, act_dim, hidden_sizes, activation, cnn=None):
        super().__init__()
        self.cnn = cnn
        self.logits_net = mlp([obs_dim] + list(hidden_sizes) + [act_dim], activation)

    def _distribution(self, obs):
        if self.cnn:
            obs = self.cnn(obs)
        logits = self.logits_net(obs)
        return Categorical(logits=logits)

    def _log_prob_from_distribution(self, pi, act):
        return pi.log_prob(act)


class GaussianActor(Actor):
    def __init__(self, obs_dim, act_dim, hidden_sizes, activation, cnn=None):
        super().__init__()
        self.cnn = cnn
        log_std = -0.5 * np.ones(act_dim, dtype=np.float32)
        self.log_std = torch.nn.Parameter(torch.as_tensor(log_std))
        self.mu_net = mlp([obs_dim] + list(hidden_sizes) + [act_dim], activation)

    def _distribution(self, obs):
        if self.cnn:
            obs = self.cnn(obs)
        mu = self.mu_net(obs)
        std = torch.exp(self.log_std)
        return Normal(mu, std)

    def _log_prob_from_distribution(self, pi, act):
        return pi.log_prob(act).sum(axis=-1)


class Critic(nn.Module):
    def __init__(self, obs_dim, hidden_sizes, activation, cnn=None):
        super().__init__()
        self.cnn = cnn
        self.v_net = mlp([obs_dim] + list(hidden_sizes) + [1], activation)

    def forward(self, obs):
        if self.cnn:
            obs = self.cnn(obs)
        return torch.squeeze(self.v_net(obs), -1)


class ActorCritic(nn.Module):
    def __init__(
        self,
        observation_dim,
        action_space,
        hidden_sizes=(64, 64),
        activation=nn.ReLU,
        cnn_enable=False,
        frames=3
    ):
        super().__init__()

        if cnn_enable:
            cnnet = cnn(frames, activation)
        else:
            cnnet = None

        obs_dim = observation_dim

        if isinstance(action_space, Box):
            self.pi = GaussianActor(obs_dim, action_space.shape[0], hidden_sizes, activation, cnnet)
        elif isinstance(action_space, Discrete):
            self.pi = CategoricalActor(obs_dim, action_space.n, hidden_sizes, activation, cnnet)
        else:
            raise NotImplementedError(f"Unsupported action space: {action_space}")

        self.v = Critic(obs_dim, hidden_sizes, activation, cnnet)

    def step(self, obs):
        with torch.no_grad():
            pi = self.pi._distribution(obs)
            a = pi.sample()
            logp_a = self.pi._log_prob_from_distribution(pi, a)
            v = self.v(obs)
        return a.cpu().numpy(), v.cpu().numpy(), logp_a.cpu().numpy()


class FramePreStacking(gym.Wrapper):
    """
    与训练代码保持一致的像素输入预处理。
    仅当训练时使用 --from_pixel 时，部署时也需要加 --from_pixel。
    """
    def __init__(self, env, n_frames=3):
        super().__init__(env)
        self.n_frames = n_frames
        self.frames = np.zeros((n_frames, *(env.observation_space.shape[:-1])))
        self.observation_space = gym.spaces.Box(
            low=0,
            high=1,
            shape=(n_frames, *(env.observation_space.shape[:-1])),
            dtype=np.float64
        )
        self.diffM = np.array([
            [1, 0, 0],
            [-2, -1, 0],
            [1, 1, 1]
        ])

    def step(self, action):
        obs, reward, terminated, truncated, info = self.env.step(action)
        self.frames = np.roll(self.frames, shift=-1, axis=0)
        self.frames[-1] = obs[:, :, 0] / 255.0
        return self.preprocess(self.frames), reward, terminated, truncated, info

    def reset(self, **kwargs):
        obs, info = self.env.reset(**kwargs)
        self.frames = np.zeros((self.n_frames, *(self.env.observation_space.shape[:-1])))
        self.frames[-1] = obs[:, :, 0] / 255.0
        return self.preprocess(self.frames), info

    def preprocess(self, I):
        I = I.transpose(1, 2, 0)
        I = np.dot(I, self.diffM)
        I = I.transpose(2, 0, 1)
        return I


def load_model(model_path, device):
    """
    加载训练保存的完整 ActorCritic 对象。
    新版 PyTorch 若遇到 weights_only 相关报错，显式 weights_only=False。
    """
    try:
        model = torch.load(model_path, map_location=device, weights_only=False)
    except TypeError:
        model = torch.load(model_path, map_location=device)

    model.to(device)
    model.eval()
    return model


def select_action(ac, obs, action_space, device, stochastic=False):
    """
    stochastic=False：部署时使用确定性动作
        - 离散动作：取概率最大的动作 argmax
        - 连续动作：取 Normal 分布均值 mean
    stochastic=True：按策略分布随机采样，便于观察策略多样性
    """
    obs_tensor = torch.as_tensor(obs, dtype=torch.float32, device=device).unsqueeze(0)

    with torch.no_grad():
        pi = ac.pi._distribution(obs_tensor)

        if stochastic:
            action_tensor = pi.sample()
        else:
            if isinstance(pi, Categorical):
                action_tensor = torch.argmax(pi.probs, dim=-1)
            elif isinstance(pi, Normal):
                action_tensor = pi.mean
            else:
                raise NotImplementedError(f"Unsupported distribution: {type(pi)}")

    action = action_tensor.cpu().numpy().squeeze()

    if isinstance(action_space, Discrete):
        return int(action)

    if isinstance(action_space, Box):
        action = np.asarray(action, dtype=np.float32)
        return np.clip(action, action_space.low, action_space.high)

    raise NotImplementedError(f"Unsupported action space: {action_space}")


def run_policy(
    env_name,
    model_path,
    device,
    episodes=5,
    max_ep_len=10000,
    from_pixel=False,
    stochastic=False,
    render_mode="human",
    sleep=0.0
):
    env = gym.make(env_name, render_mode=render_mode)

    if from_pixel:
        env = FramePreStacking(env)

    ac = load_model(model_path, device)

    for ep in range(1, episodes + 1):
        obs, info = env.reset()
        ep_ret = 0.0
        ep_len = 0

        while True:
            action = select_action(
                ac=ac,
                obs=obs,
                action_space=env.action_space,
                device=device,
                stochastic=stochastic
            )

            obs, reward, terminated, truncated, info = env.step(action)

            ep_ret += float(reward)
            ep_len += 1

            if sleep > 0:
                time.sleep(sleep)

            timeout = ep_len >= max_ep_len
            done = terminated or truncated or timeout

            if done:
                print(f"Episode {ep}: return={ep_ret:.2f}, length={ep_len}")
                break

    env.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--env", type=str, default="CartPole-v1", help="Gymnasium 环境名")
    parser.add_argument("--model", type=str, default="ppo_model.pt", help="训练保存的 pt 文件路径")
    parser.add_argument("--device", type=str, default="cpu", help="cpu 或 cuda")
    parser.add_argument("--episodes", type=int, default=5, help="部署测试回合数")
    parser.add_argument("--max_ep_len", type=int, default=10000, help="单回合最大步数")
    parser.add_argument("--from_pixel", action="store_true", help="训练时如果用了像素输入，部署时也必须开启")
    parser.add_argument("--stochastic", action="store_true", help="使用随机采样动作；默认使用确定性动作")
    parser.add_argument("--no_render", action="store_true", help="关闭窗口渲染")
    parser.add_argument("--sleep", type=float, default=0.0, help="每步暂停秒数，便于慢速观察")
    args = parser.parse_args()

    device = torch.device(args.device)
    render_mode = None if args.no_render else "human"

    run_policy(
        env_name=args.env,
        model_path=args.model,
        device=device,
        episodes=args.episodes,
        max_ep_len=args.max_ep_len,
        from_pixel=args.from_pixel,
        stochastic=args.stochastic,
        render_mode=render_mode,
        sleep=args.sleep
    )


if __name__ == "__main__":
    main()

#运行方式：
#python deploy_ppo_render.py --env CartPole-v1 --model ppo_model.pt --episodes 5 --device cpu
#如果你训练时用了像素输入：
#python deploy_ppo_render.py --env Pong-v0 --model ppo_model.pt --from_pixel --episodes 3
#如果想看随机策略表现
#python deploy_ppo_render.py --env CartPole-v1 --model ppo_model.pt --episodes 5 --stochastic