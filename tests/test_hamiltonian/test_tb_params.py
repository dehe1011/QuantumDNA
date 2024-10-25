import shutil
import pytest
import os

from qDNA import DATA_DIR
from qDNA.hamiltonian import (
    save_tb_params,
    load_tb_params,
    wrap_load_tb_params,
    wrap_save_tb_params,
)


@pytest.mark.parametrize(
    "tb_params, metadata, folder",
    [
        (
            {"t_AB": 5, "t_AC": 3, "t_BC": -2},
            {"source": "author2024", "particle": "particle", "tb_model_name": "model"},
            "delete_this_folder",
        )
    ],
)
def test_tb_params(tb_params, metadata, folder):
    save_tb_params(tb_params, metadata, folder)
    loaded_tb_params = load_tb_params(metadata, folder, load_metadata=False)
    shutil.rmtree(folder)  # Clean up test folder
    assert tb_params == loaded_tb_params


@pytest.mark.parametrize(
    "tb_params, source, particle, tb_model_name",
    [({"t_AB": 5, "t_AC": 3, "t_BC": -2}, "author2024", "particle", "model")],
)
def test_wrap_tb_params(tb_params, source, particle, tb_model_name):
    wrap_save_tb_params(tb_params, source, particle, tb_model_name)
    loaded_tb_params = wrap_load_tb_params(
        source, particle, tb_model_name, load_metadata=False
    )

    directory = os.path.join(DATA_DIR, "raw", "tb_params")
    filename = source + "_" + particle + "_" + tb_model_name + ".json"
    os.remove(os.path.join(directory, filename))
    assert tb_params == loaded_tb_params
