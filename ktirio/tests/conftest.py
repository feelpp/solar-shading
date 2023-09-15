import sys
import pytest
import feelpp

class InitFeelpp:
    def __init__(self,config):
        try:
            sys.argv=['test_sm']
            self.feelpp_env = feelpp.Environment(sys.argv,config=config)
        except Exception:
            return 


@pytest.fixture(scope="session")
def init_feelpp():
    return InitFeelpp(feelpp.globalRepository("shading-tests")).feelpp_env

@pytest.fixture(scope="session")
def init_feelpp_config_local():
    return InitFeelpp(feelpp.localRepository("feelppdb")).feelpp_env

