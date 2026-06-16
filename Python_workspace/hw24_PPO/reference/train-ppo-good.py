import numpy as np
#import gym
import gymnasium as gym
import time
from gymnasium.spaces import Box, Discrete
#from gym.spaces import Box, Discrete
import torch
import torch.nn as nn
from torch.optim import Adam
from torch.distributions.normal import Normal
from torch.distributions.categorical import Categorical

# -*- coding: utf-8 -*-

def combined_shape(length, shape=None):
    if shape is None:
        return (length,)
    return (length, shape) if np.isscalar(shape) else (length, *shape)

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
    for j in range(len(sizes)-1):
        act = activation if j < len(sizes)-2 else output_activation
        layers += [nn.Linear(sizes[j], sizes[j+1]), act()]
    return nn.Sequential(*layers)

def count_vars(module):
    return sum([np.prod(p.shape) for p in module.parameters()])

def discount_cumsum(x, discount):
    res = []
    for a in x[::-1]: # traverse from the end
        if len(res) == 0:
            res.append(a)
        else:
            res.append(a + discount * res[-1])
    return res[::-1]


class Actor(nn.Module):
    def _distribution(self, obs):  #基类
        raise NotImplementedError

    def _log_prob_from_distribution(self, pi, act):
        raise NotImplementedError

    def forward(self, obs, act=None):
        # Produce action distributions for given observations, and 
        # optionally compute the log likelihood of given actions under
        # those distributions.
        pi = self._distribution(obs)
        logp_a = None
        if act is not None:
            logp_a = self._log_prob_from_distribution(pi, act)
        return pi, logp_a 


class CategoricalActor(Actor):  #离散动作空间,指定继承Actor基类
    def __init__(self, obs_dim, act_dim, hidden_sizes, activation, cnn=None):
        super().__init__()
        self.cnn = cnn
        self.logits_net = mlp([obs_dim] + list(hidden_sizes) + [act_dim], activation)
        #log(p / (1 - p))
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
        self.log_std = torch.nn.Parameter(torch.as_tensor(log_std))  # learnable parameter,not a network
        #torch.nn.Parameter will be registered as model parameter
        #self.log_std 不是一个网络，而是一个可学习参数
        #它是标量或向量，用于定义动作分布的标准差
        #通过 torch.nn.Parameter 使其成为模型参数，可以在训练中优化
        self.mu_net = mlp([obs_dim] + list(hidden_sizes) + [act_dim], activation)
        # input: observation state, output: mean of act distribution
    def _distribution(self, obs):
        if self.cnn:
            obs = self.cnn(obs)
        mu = self.mu_net(obs)
        std = torch.exp(self.log_std)
        return Normal(mu, std)

    def _log_prob_from_distribution(self, pi, act):
        return pi.log_prob(act).sum(axis=-1)    # Last axis sum needed for Torch Normal distribution
        # log_prob: log probability of action under the distribution
        #这是策略梯度方法中计算重要性采样比率的关键组成部分
        #sum(axis=-1) 部分：对动作的所有维度求和，得到整个动作的联合对数概率
        #这是因为多维正态分布的联合概率密度等于各维度概率密度的乘积
        #在对数空间中，乘积变为求和
class Critic(nn.Module):
    def __init__(self, obs_dim, hidden_sizes, activation, cnn=None):
        super().__init__()
        self.cnn = cnn
        self.v_net = mlp([obs_dim] + list(hidden_sizes) + [1], activation)

    def forward(self, obs):
        if self.cnn:
            obs = self.cnn(obs)
        return torch.squeeze(self.v_net(obs), -1)
        # Critical to ensure v has right shape.


class ActorCritic(nn.Module):
    def __init__(self, observation_dim, action_space, 
                 hidden_sizes=(64,64), activation=nn.ReLU, cnn_enable=False, frames=3):
        super().__init__()

        if cnn_enable:
            cnnet = cnn(frames, activation)
        else:
            cnnet = None
        obs_dim = observation_dim

        # policy builder depends on action space
        if isinstance(action_space, Box):
            self.pi = GaussianActor(obs_dim, action_space.shape[0], hidden_sizes, activation, cnnet)
            #检查动作空间是否为 Box 类型，这是 Gym 中表示连续动作空间的类。例如，动作可以是任意实数值
        elif isinstance(action_space, Discrete):
            self.pi = CategoricalActor(obs_dim, action_space.n, hidden_sizes, activation, cnnet)
            #检查动作空间是否为 Discrete 类型，这是 Gym 中表示离散动作空间的类。
            # 例如，动作只能从有限的选项中选择，如{0, 1, 2, 3}
        # build value function
        self.v  = Critic(obs_dim, hidden_sizes, activation, cnnet)

    def step(self, obs):
        #接受观测值 obs 作为输入，用于在环境中执行一步操作
        with torch.no_grad(): #不需要计算梯度，因为这是在环境中执行动作，而不是训练模型
            pi = self.pi._distribution(obs)
            a = pi.sample()
            logp_a = self.pi._log_prob_from_distribution(pi, a)
            v = self.v(obs)
            #调用价值网络（critic）计算当前观测状态的价值估计 v
        return a.cpu().numpy(), v.cpu().numpy(), logp_a.cpu().numpy()
        #返回动作 a，价值估计 v，动作概率对数 logp_a

class FramePreStacking(gym.Wrapper):  #帧预处理,用于处理基于像素的强化学习环境
    #继承自 gym.Wrapper，是一个 Gym 环境包装器基类
    def __init__(self, env, n_frames=3): #默认堆叠3帧图像
        super().__init__(env)
        self.n_frames = n_frames
        self.frames = np.zeros((n_frames, *(env.observation_space.shape[:-1])))
        #存储最近的n帧图像数据
        self.observation_space = gym.spaces.Box(low=0, high=1, shape=(n_frames, *(env.observation_space.shape[:-1])), dtype=np.float64)
        #定义新的观测空间为堆叠后的桢，形状为(n_frames, height, width)，值范围在0到1之间
        self.diffM = np.array([[1 ,0 ,0],  #位置：当前帧
                               [-2,-1,0],  # 速度：当前帧-前一帧
                               [1 ,1 ,1]]) # 加速度：速度变化

    def step(self, action):
        obs, reward, done, truncated, info = self.env.step(action) 
        #调用原始环境的 step 方法获取当前帧
        self.frames = np.roll(self.frames, shift=-1, axis=0)
        #将旧帧向前移动（最旧的帧被移除）
        self.frames[-1] = obs[:,:,0]/255.
        #将新观测帧添加到最后一层，并归一化到[0,1]范围
        return self.preprocess(self.frames), reward, done, truncated, info
        #返回堆叠后的桢作为新的观测值， reward，done，truncated，info
    def reset(self):
        obs, _ = self.env.reset()
        self.frames = np.zeros((self.n_frames, *(self.env.observation_space.shape[:-1])))
        self.frames[-1] = obs[:,:,0]/255.
        return self.preprocess(self.frames), _
        #返回堆叠后的桢作为新的观测值， reward，done，truncated，info
        #_: 从原始环境 reset() 方法中获得的额外信息，原样传递
    def preprocess(self, I):
        # compute position, velocity, acceleration
        I = I.transpose(1,2,0) ## 改变维度顺序
        #原来的第1维变为新的第0维，原来的第2维变为新的第1维，原来的第0维变为新的第2维，
        I = np.dot(I, self.diffM)
        I = I.transpose(2,0,1) ## 恢复维度顺序
        return I


class PPOBuffer:
    """
    A buffer for storing trajectories experienced by a PPO agent interacting
    with the environment, and using Generalized Advantage Estimation (GAE-Lambda)
    for calculating the advantages of state-action pairs.
    """
    def __init__(self, obs_dim, act_dim, size, gamma=0.99, lam=0.95, device='cpu'):
        self.obs_buf = np.zeros(combined_shape(size, obs_dim), dtype=np.float32)
        self.act_buf = np.zeros(combined_shape(size, act_dim), dtype=np.float32)
        self.adv_buf = np.zeros(size, dtype=np.float32)  # 存储每个时间步的优势advantage
        self.rew_buf = np.zeros(size, dtype=np.float32)  # 存储每个时间步的奖励reward
        self.ret_buf = np.zeros(size, dtype=np.float32)  # 存储每个状态的回报值（Return/Discounted Reward-to-go）
            #R_t = r_t + γ*r_{t+1} + γ^2*r_{t+2} + ... 
            #ret_buf用于训练价值网络，使其能准确估计 V(s)
            #adv_buf用于训练策略网络，指导策略改进方向
        self.val_buf = np.zeros(size, dtype=np.float32) #状态价值估计
        self.logp_buf = np.zeros(size, dtype=np.float32) #动作的对数概率
        self.gamma, self.lam = gamma, lam
        self.ptr, self.path_start_idx, self.max_size = 0, 0, size
        self.device = device

    def store(self, obs, act, rew, val, logp):
        """
        Append one timestep of agent-environment interaction to the buffer.
        """
        assert self.ptr < self.max_size     # buffer has to have room so you can store
        self.obs_buf[self.ptr] = obs
        self.act_buf[self.ptr] = act
        self.rew_buf[self.ptr] = rew
        self.val_buf[self.ptr] = val
        self.logp_buf[self.ptr] = logp
        self.ptr += 1

    def finish_path(self, last_val=0):
        """
        Call this at the end of a trajectory, or when one gets cut off
        by an epoch ending. This looks back in the buffer to where the
        trajectory started, and uses rewards and value estimates from
        the whole trajectory to compute advantage estimates with GAE-Lambda,
        as well as compute the rewards-to-go for each state, to use as
        the targets for the value function.

        The "last_val" argument should be 0 if the trajectory ended
        because the agent reached a terminal state (died), and otherwise
        should be V(s_T), the value function estimated for the last state.
        This allows us to bootstrap the reward-to-go calculation to account
        for timesteps beyond the arbitrary episode horizon (or epoch cutoff).
        """

        path_slice = slice(self.path_start_idx, self.ptr)
        rews = np.append(self.rew_buf[path_slice], last_val)
        vals = np.append(self.val_buf[path_slice], last_val)
        
        # the next two lines implement GAE-Lambda advantage calculation
        deltas = rews[:-1] + self.gamma * vals[1:] - vals[:-1]
        self.adv_buf[path_slice] = discount_cumsum(deltas, self.gamma * self.lam)
        
        # the next line computes rewards-to-go, to be targets for the value function
        self.ret_buf[path_slice] = discount_cumsum(rews, self.gamma)[:-1]
        # another way to compute return: R = A + V (better to be a moving average of V)
        # self.ret_buf[path_slice] = self.adv_buf[path_slice] + self.val_buf[path_slice]
        
        self.path_start_idx = self.ptr

    def get(self):  #在每个epoch结束时调用，获取所有存储的数据
        """
        Call this at the end of an epoch to get all of the data from
        the buffer, with advantages appropriately normalized (shifted to have
        mean zero and std one). Also, resets some pointers in the buffer.
        """
        assert self.ptr == self.max_size    # buffer has to be full before you can get
        self.ptr, self.path_start_idx = 0, 0
        # the next two lines implement the advantage normalization trick
        adv_mean, adv_std = self.adv_buf.mean(), self.adv_buf.std()
        self.adv_buf = (self.adv_buf - adv_mean) / (adv_std + 1e-8)
        data = dict(obs=self.obs_buf, act=self.act_buf, ret=self.ret_buf,
                    adv=self.adv_buf, logp=self.logp_buf)
        return {k: torch.as_tensor(v, dtype=torch.float32).to(self.device) for k,v in data.items()}


def ppo(env_fn, actor_critic=ActorCritic, ac_kwargs=dict(), seed=0, 
        steps_per_epoch=6000, epochs=50, gamma=0.99, clip_ratio=0.2, pi_lr=3e-3,
        vf_lr=1e-4, train_pi_iters=10, train_v_iters=80, lam=0.97, max_ep_len=6000,
        target_kl=10, save_freq=10, load_from=None, device='cpu'):
    #pi_lr=3e-4，target_kl=0.1，train_v_iters=10，vf_lr=3e-4
    print(locals())

    #with open(__file__, 'rb') as f:
    #    print(f.read())
    #打开并读取当前脚本文件（即 ppo.py 文件）的全部内容，并打印出来
    # Random seed
    seed = 42
    torch.manual_seed(seed)
    np.random.seed(seed)

    # Instantiate environment
    env = env_fn()
    act_dim = env.action_space.shape

    if args.from_pixel:
        cnn_enable = True
        env = FramePreStacking(env)
        obs_dim = (13824) # a hack number for 210x160 image and our cnn.
    else:
        cnn_enable = False
        obs_dim = env.observation_space.shape[0]

    # Create actor-critic module
    if not load_from:
        ac = actor_critic(obs_dim, env.action_space, cnn_enable=cnn_enable, **ac_kwargs).to(device)
    else:
        ac = torch.load('ppo_model.pt', map_location=device)

    # Count variables
    var_counts = tuple(count_vars(module) for module in [ac.pi, ac.v])
    print('\nNumber of parameters: \t pi: %d, \t v: %d\n'%var_counts)

    # Set up experience buffer
    buf = PPOBuffer(env.observation_space.shape, act_dim, steps_per_epoch, gamma, lam, device)

    # Set up function for computing PPO policy loss
    def compute_loss_pi(data):
        obs, act, adv, logp_old = data['obs'], data['act'], data['adv'], data['logp']

        # Policy loss
        pi, logp = ac.pi(obs, act)
        ratio = torch.exp(logp - logp_old)
        clip_adv = torch.clamp(ratio, 1-clip_ratio, 1+clip_ratio) * adv
        loss_pi = -(torch.min(ratio * adv, clip_adv)).mean()
        # Useful extra info
        approx_kl = (logp_old - logp).mean().item()  #近似KL散度，.item(): 将张量转换为Python标量值
        ent = pi.entropy().mean().item() # you can also add entropy bonus to loss_pi to encourage exploration.
        clipped = ratio.gt(1+clip_ratio) | ratio.lt(1-clip_ratio)
        clipfrac = torch.as_tensor(clipped, dtype=torch.float32).mean().item()
        pi_info = dict(kl=approx_kl, ent=ent, cf=clipfrac)

        return loss_pi, pi_info

    # Set up function for computing value loss, you can also use clipped version of loss_v.
    def compute_loss_v(data):
        obs, ret = data['obs'], data['ret']
        return ((ac.v(obs) - ret)**2).mean()

    # Set up optimizers for policy and value function
    pi_optimizer = Adam(ac.pi.parameters(), lr=pi_lr)
    vf_optimizer = Adam(ac.v.parameters(), lr=vf_lr)


    def update():
        data = buf.get() # here we use a single batch of data, you can also use multiple mini-batches of data to update.
        # Train policy with multiple steps of gradient descent
        for i in range(train_pi_iters):
            pi_optimizer.zero_grad()
            loss_pi, pi_info = compute_loss_pi(data)
            kl = pi_info['kl']
            if kl > 1.5 * target_kl: # early stopping 限制 KL 散度，防止策略更新过大导致性能下降
                print('Early stopping at step %d due to reaching max kl.'%i)
                break
            loss_pi.backward()
            pi_optimizer.step()


        # Value function learning
        for i in range(train_v_iters):
            vf_optimizer.zero_grad()
            loss_v = compute_loss_v(data)
            loss_v.backward()
            vf_optimizer.step()

        # Log changes from update
        kl, ent, cf = pi_info['kl'], pi_info['ent'], pi_info['cf']
        print(dict(LossPi=loss_pi.item(), LossV=loss_v.item(),
                     KL=kl, Entropy=ent, ClipFrac=cf))

    # Prepare for interaction with environment
    start_time = time.time()
    o, _ = env.reset()
    ep_ret, ep_len = 0, 0

    # Main loop: collect experience in env and update/log each epoch
    for epoch in range(epochs):
        print('Epoch: {}'.format(epoch))
        for t in range(steps_per_epoch):
            a, v, logp = ac.step(torch.as_tensor(o, dtype=torch.float32).unsqueeze(0).to(device))
            a = a.squeeze(0)

            next_o, r, d, _, _ = env.step(a)
            ep_ret += r
            ep_len += 1

            buf.store(o, a, r, v, logp)
            o = next_o

            timeout = ep_len == max_ep_len
            terminal = d or timeout
            epoch_ended = t==steps_per_epoch-1
            pong_ended = (r != 0 and 'Pong' in args.env)

            if terminal or epoch_ended or pong_ended:
                # if epoch_ended and not(terminal):
                #     print('Warning: trajectory cut off by epoch at %d steps.'%ep_len, flush=True)
                # if trajectory didn't reach terminal state, bootstrap value target
                if timeout or epoch_ended:
                    _, v, _ = ac.step(torch.as_tensor(o, dtype=torch.float32).unsqueeze(0).to(device))
                else:
                    v = 0
                buf.finish_path(v)
                # if terminal:
                #     print(dict(EpRet=ep_ret, EpLen=ep_len))
                if terminal or epoch_ended:
                    print(dict(EpRet=ep_ret, EpLen=ep_len))
                    o, _ = env.reset()
                    ep_ret, ep_len = 0, 0


        # Save model
        if (epoch % save_freq == 0) or (epoch == epochs-1):
            torch.save(ac, 'ppo_model.pt')

        # Perform PPO update!
        update()
    print('Total time: {:.2f}s'.format(time.time() - start_time))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument('--env', type=str, default='Pong-v0')
    parser.add_argument('--env', type=str, default='CartPole-v1')
    parser.add_argument('--hid', type=int, default=256)
    parser.add_argument('--l', type=int, default=1)  # number of layers
    parser.add_argument('--kl', type=float, default=0.01)  # KL
    parser.add_argument('--gamma', type=float, default=0.99)  # discount factor
    parser.add_argument('--seed', '-s', type=int, default=0)  # random seed
    #parser.add_argument('--steps', type=int, default=6000) 
    parser.add_argument('--steps', type=int, default=500)   
    parser.add_argument('--epochs', type=int, default=10000)   
    parser.add_argument('--exp_name', type=str, default='ppo')    
    #parser.add_argument('--device', type=str, default='cuda')
    parser.add_argument('--device', type=str, default='cpu')
    parser.add_argument('--from_pixel', action="store_true") #store_true,当命令行中出现该参数时，其值为 True；当不出现时，其值为 False。
    parser.add_argument('--render', action="store_true",default=True) #default=True,渲染
    #parser.add_argument('--render', action="store_true")
    parser.add_argument('--load_from', action="store_true",default=True) #False
    #parser.add_argument('--load_from', action="store_true")
    args = parser.parse_args()
    '''
    if args.env.startswith('Pong'):
        args.from_pixel = True
    '''
    ppo(lambda : gym.make(args.env, render_mode='human' if args.render else None), actor_critic=ActorCritic,
        ac_kwargs=dict(hidden_sizes=[args.hid]*args.l), target_kl=args.kl, gamma=args.gamma, 
        seed=args.seed, steps_per_epoch=args.steps, epochs=args.epochs, load_from=args.load_from, device=torch.device(args.device))

    '''
    这是整个PPO算法的入口点,它将启动完整的训练流程，包括：
    环境交互、经验收集、策略更新、价值函数更新、模型保存
    使用lambda表达式创建Gym环境实例
    args.env:从命令行参数获取环境名称（默认是'CartPole-v1')
    render_mode:根据args.render参数决定是否渲染环境
    如果args.render为True,则render_mode='human'（可视化显示）,否则为None(无图形界面运行）
    ac_kwargs:传递给ActorCritic构造函数的参数字典 (kwargs=keyword arguments)
    target_kl=args.kl,      # KL散度限制
    gamma=args.gamma,       # 折扣因子
    seed=args.seed,         # 随机种子
    steps_per_epoch=args.steps,  # 每个epoch的步数
    epochs=args.epochs,     # 训练epoch数
    load_from=args.load_from,    # 是否从已有模型加载
    '''