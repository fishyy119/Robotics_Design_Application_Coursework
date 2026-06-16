from typing import TypeAlias

import numpy as np
from jaxtyping import Float, Int
from torch import Tensor

ObsArray: TypeAlias = Float[np.ndarray, "obs_dim"]
ObsBatchArray: TypeAlias = Float[np.ndarray, "steps obs_dim"]
ScalarArray: TypeAlias = Float[np.ndarray, "steps"]
ActionArray: TypeAlias = Int[np.ndarray, "steps"]

ObsBatchTensor: TypeAlias = Float[Tensor, "batch obs_dim"]
ActionTensor: TypeAlias = Int[Tensor, "batch"]
ValueTensor: TypeAlias = Float[Tensor, "batch"]
LogitsTensor: TypeAlias = Float[Tensor, "batch act_dim"]
