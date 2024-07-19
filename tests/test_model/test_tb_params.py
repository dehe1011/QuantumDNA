import pytest
from qDNA.model import save_tb_params, load_tb_params

def test_tb_params():
    tb_param_dict = {'relation_AliceBob': 5, 'relation_AliceCharlie': 3, 'relation_AliceEve': -2,
            'relation_BobAlice': 5, 'relation_CharlieAlice': 3, 'relation_EveAlice': -2}
    info_dict={'author': 'Herb2024', 'subject': 'relations_between_persons'}
    save_tb_params(tb_param_dict, info_dict, directory = 'data/raw/test_params', notes = 'The parameters describe relations between persons.')
    loaded_tb_param_dict = load_tb_params(info_dict, directory = 'data/raw/test_params', load_metadata = False)
    assert tb_param_dict == loaded_tb_param_dict

if __name__ == "__main__":
    pytest.main()
