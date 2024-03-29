from io import StringIO
from typing import Union, Any, Type

import numpy as np

from .xlsx import dict_to_xlsx
from .. import dists


class InputParser:
    """Converts """

    def __init__(self, dist_params: dict, n: int):
        assert isinstance(dist_params, dict)
        assert isinstance(n, int)

        self.__n = n
        self.__in_raw = dist_params
        self.__in = InputParser.unflatten_dict(dist_params)

    def to_dict(self):
        n = self.__n
        dist_params = self.__in
        dict_out = dict()

        for k, v in dist_params.items():
            if isinstance(v, float) or isinstance(v, int) or isinstance(v, float):
                dict_out[k] = np.full((n,), v, dtype=float)
            elif isinstance(v, str):
                dict_out[k] = np.full(
                    (n,), v, dtype=np.dtype("U{:d}".format(len(v)))
                )
            elif isinstance(v, np.ndarray) or isinstance(v, list):
                dict_out[k] = list(np.full((n, len(v)), v, dtype=float))
            elif isinstance(v, dict):
                if "dist" in v:
                    try:
                        dict_out[k] = InputParser._sampling(v, n)
                    except KeyError:
                        raise KeyError(f"Missing parameters in input variable {k}.")
                elif "ramp" in v:
                    s_ = StringIO(v["ramp"])
                    d_ = np.loadtxt(s_, delimiter=',')
                    t_ = d_[:, 0]
                    v_ = d_[:, 1]
                    if all(v_ == v_[0]):
                        f_interp = v_[0]
                    else:
                        def f_interp(x):
                            return np.interp(x, t_, v_)
                    dict_out[k] = np.full((n,), f_interp)
                else:
                    raise ValueError(f"Unknown input data type for {k}. {v}.")
            elif v is None:
                dict_out[k] = np.full((n,), np.nan, dtype=float)
            else:
                raise TypeError(f"Unknown input data type for {k}.")

        dict_out["index"] = np.arange(0, n, 1)
        return dict_out

    def to_xlsx(self, fp: str):
        dict_to_xlsx({i: InputParser.flatten_dict(v) for i, v in self.to_dict().items()}, fp)

    @staticmethod
    def unflatten_dict(dict_in: dict) -> dict:
        """Invert flatten_dict.

        :param dict_in:
        :return dict_out:
        """
        dict_out = dict()

        for k, v in dict_in.items():
            InputParser.__unflatten_dict(k, v, dict_out)

        return dict_out

    @staticmethod
    def __unflatten_dict(k: str, v: Any, dict_out: dict):
        if ":" in k:
            k1, *k2 = k.split(':')
            if k1 not in dict_out:
                dict_out[k1] = dict()
            InputParser.__unflatten_dict(':'.join(k2), v, dict_out[k1])
        else:
            dict_out[k] = v

    @staticmethod
    def flatten_dict(dict_in: dict) -> dict:
        dict_out = dict()
        InputParser.__flatten_dict(dict_in, dict_out)
        return dict_out

    @staticmethod
    def __flatten_dict(dict_in: dict, dict_out: dict, history: str = None):
        """Converts two levels dict to single level dict. Example input and output see _test_dict_flatten.
        >>> dict_in = {
        >>>             'a': 1,
        >>>             'b': {'b1': 21, 'b2': 22},
        >>>             'c': {'c1': 31, 'c2': 32, 'c3': 33}
        >>>         }
        >>> output = {
        >>>             'a': 1,
        >>>             'b:b1': 21,
        >>>             'b:b2': 22,
        >>>             'c:c1': 31,
        >>>             'c:c2': 32,
        >>>             'c:c3': 33,
        >>>         }
        >>> assert InputParser.flatten_dict(dict_in) == output  # True

        :param dict_in:     Any two levels (or less) dict.
        :return dict_out:   Single level dict.
        """
        for k, v in dict_in.items():
            if isinstance(v, dict):
                InputParser.__flatten_dict(v, dict_out=dict_out, history=k if history is None else f'{history}:{k}')
            else:
                dict_out[f'{k}' if history is None else f'{history}:{k}'] = v

    @staticmethod
    def _sampling(dist_params: dict, num_samples: int, randomise: bool = True) -> Union[float, np.ndarray]:
        """A reimplementation of _sampling_scipy but without scipy"""
        dist_name = ''.join(dist_params.pop('dist').replace('_', ' ').strip().title().split())
        if dist_name == 'Norm':
            dist_name = 'Normal'
        elif dist_name == 'GumbelR':
            dist_name = 'Gumbel'
        elif dist_name == 'Uniform':
            if (
                    'lbound' in dist_params and 'ubound' in dist_params and
                    'mean' not in dist_params and 'sd' not in dist_params
            ):
                a = dist_params.pop('lbound')
                b = dist_params.pop('ubound')
                mean = (a + b) / 2
                sd = (b - a) / (2 * np.sqrt(3))
                dist_params['mean'] = mean
                dist_params['sd'] = sd
        elif dist_name == 'LognormMod':
            dist_name = 'LognormalMod'
        elif dist_name == 'Lognorm':
            dist_name = 'Lognormal'
        elif dist_name == 'Constant':
            if 'ubound' in dist_params and 'lbound' in dist_params:
                dist_params['value'] = (dist_params.pop('lbound') + dist_params.pop('ubound')) / 2.

        dist_cls: Type[dists.DistFunc] = getattr(dists, dist_name)
        dist_obj: dists.DistFunc = dist_cls(**dist_params)
        lim_1 = None if 'lbound' not in dist_params else dist_params['lbound']
        lim_2 = None if 'ubound' not in dist_params else dist_params['ubound']
        return dist_obj.sampling(n=num_samples, lim_1=lim_1, lim_2=lim_2, shuffle=randomise)
