import reframe as rfm
import os
import reframe.utility.sanity as sn

@rfm.simple_test
class RNGBenchmark(rfm.RegressionTest):
    # Define the variables
    cpp_file_name = parameter(['uniform_int'])
    flagsArr = parameter([
        "-O1", "-O2", "-O3", "-Ofast",
        "-O1 -march=native", "-O2 -march=native", "-O3 -march=native", "-Ofast -march=native",
        "-O1 -march=native -funroll-loops", "-O2 -march=native -funroll-loops", "-O3 -march=native -funroll-loops", "-Ofast -march=native -funroll-loops",
        "-O1 -march=native -fdisable-tree-cunrolli", "-O2 -march=native -fdisable-tree-cunrolli", "-O3 -march=native -fdisable-tree-cunrolli", "-Ofast -march=native -fdisable-tree-cunrolli"
    ])

    # Set the required programming environments and systems
    valid_prog_environs = ['*']
    valid_systems = ['*']

    # Set the build system
    build_system = 'SingleSource'
    sourcesdir = 'src'

    def __init__(self):
        self.sourcepath = f'{self.cpp_file_name}.cpp'
        self.build_system = 'SingleSource'
        self.build_system.cxx = 'g++'
        self.build_system.cxxflags = [
            '-Wno-enum-compare',
            '-I' + os.path.abspath('../extlibs/EigenRand'),
            '-I' + os.path.abspath('../extlibs/Xoshiro-cpp'),
            '-I' + os.path.abspath('../extlibs/pcg-cpp/include'),
            '-I' + os.path.abspath('../extlibs/eigen')
        ]

    # Define the pre-run commands
    @run_before('run')
    def set_prerun_cmds(self):
        self.prerun_cmds = []
        for flag in self.flagsArr:
            build_cmd = f'g++ -DFROM_SHELL_SCRIPT -Wno-enum-compare {flag} -o {self.cpp_file_name} {self.sourcepath}'
            self.prerun_cmds.extend([
                f'echo "Build command: {build_cmd}"',
                f'echo "Executable file: {self.cpp_file_name}"',
                build_cmd,
                f'./{self.cpp_file_name} {flag}'
            ])

    # Define the sanity patterns
    @sanity_function
    def assert_output(self):
        return sn.all([
            sn.assert_found(r'Eigen::VectorXi::Random', self.stdout),
            sn.assert_found(r'std::uniform_int_distribution', self.stdout),
            sn.assert_found(r'Mersenne Twister', self.stdout),
            sn.assert_found(r'EigenRanda', self.stdout),
            sn.assert_found(r'EigenRandb', self.stdout),
            sn.assert_found(r'pcg-cpp', self.stdout),
            sn.assert_found(r'XoshiroCpp', self.stdout)
        ])

    # Define the performance patterns
    @performance_function('ms')
    def extract_timing_eigen(self):
        return sn.extractsingle(r'Eigen::VectorXi::Random,(\S+)', self.stdout, 1, float)
    
    @performance_function('ms')
    def extract_timing_std_uniform(self):
        return sn.extractsingle(r'std::uniform_int_distribution,(\S+)', self.stdout, 1, float)

    @performance_function('ms')
    def extract_timing_mersenne(self):
        return sn.extractsingle(r'Mersenne Twister,(\S+)', self.stdout, 1, float)

    @performance_function('ms')
    def extract_timing_eigenranda(self):
        return sn.extractsingle(r'EigenRanda,(\S+)', self.stdout, 1, float)

    @performance_function('ms')
    def extract_timing_eigenrandb(self):
        return sn.extractsingle(r'EigenRandb,(\S+)', self.stdout, 1, float)

    @performance_function('ms')
    def extract_timing_pcg_cpp(self):
        return sn.extractsingle(r'pcg-cpp,(\S+)', self.stdout, 1, float)

    @performance_function('ms')
    def extract_timing_xoshirocpp(self):
        return sn.extractsingle(r'XoshiroCpp,(\S+)', self.stdout, 1, float)
