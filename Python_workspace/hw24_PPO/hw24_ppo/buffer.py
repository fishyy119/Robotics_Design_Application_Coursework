from dataclasses import dataclass
from typing import cast

import numpy as np
import torch

from .types import (
    ActionArray,
    ActionTensor,
    ObsArray,
    ObsBatchArray,
    ObsBatchTensor,
    ScalarArray,
    ValueTensor,
)


@dataclass(frozen=True)
class RolloutBatch:
    obs: ObsBatchTensor
    act: ActionTensor
    adv: ValueTensor
    ret: ValueTensor
    logp: ValueTensor


class RolloutBuffer:
    def __init__(
        self,
        obs_dim: int,
        size: int,
        gamma: float,
        lam: float,
    ) -> None:
        self.obs_buf: ObsBatchArray = cast(ObsBatchArray, np.zeros((size, obs_dim), dtype=np.float32))
        self.act_buf: ActionArray = cast(ActionArray, np.zeros(size, dtype=np.int64))
        self.rew_buf: ScalarArray = cast(ScalarArray, np.zeros(size, dtype=np.float32))
        self.adv_buf: ScalarArray = cast(ScalarArray, np.zeros(size, dtype=np.float32))
        self.ret_buf: ScalarArray = cast(ScalarArray, np.zeros(size, dtype=np.float32))
        self.val_buf: ScalarArray = cast(ScalarArray, np.zeros(size, dtype=np.float32))
        self.logp_buf: ScalarArray = cast(ScalarArray, np.zeros(size, dtype=np.float32))
        self.gamma = gamma
        self.lam = lam
        self.max_size = size
        self.ptr = 0
        self.path_start_idx = 0

    def store(
        self,
        obs: ObsArray,
        action: int,
        reward: float,
        value: float,
        log_prob: float,
    ) -> None:
        if self.ptr >= self.max_size:
            raise RuntimeError("RolloutBuffer is full.")
        self.obs_buf[self.ptr] = obs
        self.act_buf[self.ptr] = int(action)
        self.rew_buf[self.ptr] = float(reward)
        self.val_buf[self.ptr] = float(value)
        self.logp_buf[self.ptr] = float(log_prob)
        self.ptr += 1

    def finish_path(self, last_value: float = 0.0) -> None:
        path_slice = slice(self.path_start_idx, self.ptr)
        rewards = np.append(self.rew_buf[path_slice], last_value).astype(np.float32)
        values = np.append(self.val_buf[path_slice], last_value).astype(np.float32)
        deltas = rewards[:-1] + self.gamma * values[1:] - values[:-1]
        self.adv_buf[path_slice] = discount_cumsum(cast(ScalarArray, deltas), self.gamma * self.lam)
        discounted_returns = discount_cumsum(cast(ScalarArray, rewards), self.gamma)
        self.ret_buf[path_slice] = discounted_returns[:-1]
        self.path_start_idx = self.ptr

    def get_tensors(self, device: torch.device) -> RolloutBatch:
        if self.ptr != self.max_size:
            raise RuntimeError("RolloutBuffer must be full before sampling.")

        self.ptr = 0
        self.path_start_idx = 0

        advantage_mean = float(np.mean(self.adv_buf))
        advantage_std = float(np.std(self.adv_buf))
        normalized_advantage = (self.adv_buf - advantage_mean) / (advantage_std + 1e-8)

        obs_tensor: ObsBatchTensor = torch.as_tensor(self.obs_buf, dtype=torch.float32, device=device)
        act_tensor: ActionTensor = torch.as_tensor(self.act_buf, dtype=torch.int64, device=device)
        adv_tensor: ValueTensor = torch.as_tensor(normalized_advantage, dtype=torch.float32, device=device)
        ret_tensor: ValueTensor = torch.as_tensor(self.ret_buf, dtype=torch.float32, device=device)
        logp_tensor: ValueTensor = torch.as_tensor(self.logp_buf, dtype=torch.float32, device=device)

        return RolloutBatch(
            obs=obs_tensor,
            act=act_tensor,
            adv=adv_tensor,
            ret=ret_tensor,
            logp=logp_tensor,
        )


def discount_cumsum(values: ScalarArray, discount: float) -> ScalarArray:
    discounted = np.zeros_like(values, dtype=np.float32)
    running_total = 0.0
    for index in range(values.shape[0] - 1, -1, -1):
        running_total = float(values[index]) + discount * running_total
        discounted[index] = running_total
    return cast(ScalarArray, discounted)
