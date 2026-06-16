from typing import Sequence, cast

import torch
from torch import nn
from torch.distributions.categorical import Categorical

from .config import ModelConfig
from .types import ActionTensor, LogitsTensor, ObsBatchTensor, ValueTensor


def get_activation(name: str) -> type[nn.Module]:
    activations: dict[str, type[nn.Module]] = {
        "relu": nn.ReLU,
        "tanh": nn.Tanh,
    }
    if name not in activations:
        raise ValueError(f"Unsupported activation: {name}.")
    return activations[name]


def build_mlp(
    layer_sizes: Sequence[int],
    activation_name: str,
    output_activation: type[nn.Module] = nn.Identity,
) -> nn.Sequential:
    layers: list[nn.Module] = []
    activation = get_activation(activation_name)
    for index in range(len(layer_sizes) - 1):
        input_dim = int(layer_sizes[index])
        output_dim = int(layer_sizes[index + 1])
        current_activation = activation if index < len(layer_sizes) - 2 else output_activation
        layers.append(nn.Linear(input_dim, output_dim))
        layers.append(current_activation())
    return nn.Sequential(*layers)


class CategoricalActor(nn.Module):
    def __init__(
        self,
        obs_dim: int,
        act_dim: int,
        hidden_sizes: Sequence[int],
        activation_name: str,
    ) -> None:
        super().__init__()
        layer_sizes = [obs_dim, *hidden_sizes, act_dim]
        self.logits_net = build_mlp(layer_sizes, activation_name)

    def forward(self, obs: ObsBatchTensor) -> LogitsTensor:
        logits = self.logits_net(obs)
        return cast(LogitsTensor, logits)

    def distribution(self, obs: ObsBatchTensor) -> Categorical:
        return Categorical(logits=self.forward(obs))


class Critic(nn.Module):
    def __init__(
        self,
        obs_dim: int,
        hidden_sizes: Sequence[int],
        activation_name: str,
    ) -> None:
        super().__init__()
        layer_sizes = [obs_dim, *hidden_sizes, 1]
        self.value_net = build_mlp(layer_sizes, activation_name)

    def forward(self, obs: ObsBatchTensor) -> ValueTensor:
        values = self.value_net(obs).squeeze(-1)
        return cast(ValueTensor, values)


class MLPActorCritic(nn.Module):
    def __init__(self, obs_dim: int, act_dim: int, config: ModelConfig) -> None:
        super().__init__()
        self.actor = CategoricalActor(
            obs_dim=obs_dim,
            act_dim=act_dim,
            hidden_sizes=config.hidden_sizes,
            activation_name=config.activation,
        )
        self.critic = Critic(
            obs_dim=obs_dim,
            hidden_sizes=config.hidden_sizes,
            activation_name=config.activation,
        )

    def act(self, obs: ObsBatchTensor) -> tuple[ActionTensor, ValueTensor, ValueTensor]:
        with torch.no_grad():
            distribution = self.actor.distribution(obs)
            actions = distribution.sample()
            log_prob = cast(ValueTensor, distribution.log_prob(actions))
            values = self.critic(obs)
        return actions, values, log_prob

    def act_deterministic(self, obs: ObsBatchTensor) -> ActionTensor:
        with torch.no_grad():
            logits = self.actor(obs)
            actions = torch.argmax(logits, dim=-1)
        return actions

    def evaluate_actions(
        self, obs: ObsBatchTensor, actions: ActionTensor
    ) -> tuple[ValueTensor, ValueTensor, ValueTensor]:
        distribution = self.actor.distribution(obs)
        log_prob = cast(ValueTensor, distribution.log_prob(actions))
        entropy = distribution.entropy()
        values = self.critic(obs)
        return log_prob, entropy, values
