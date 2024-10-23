import shutil
from qDNA.model import save_tb_params, load_tb_params


def test_tb_params():
    tb_params = {"t_AB": 5, "t_AC": 3, "t_BC": -2}
    metadata = {
        "source": "author2024",
        "particle": "particle",
        "tb_model_name": "model",
    }
    save_tb_params(tb_params, metadata, "delete_this_folder")
    loaded_tb_params = load_tb_params(
        metadata, "delete_this_folder", load_metadata=False
    )
    shutil.rmtree("delete_this_folder")
    assert tb_params == loaded_tb_params
