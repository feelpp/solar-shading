import reframe as rfm
import os
import reframe.utility.sanity as sn

flags = [
    "-O1", "-O2", "-O3", "-Ofast",
    "-O1 -march=native", "-O2 -march=native", "-O3 -march=native", "-Ofast -march=native",
    "-O1 -march=native -funroll-loops", "-O2 -march=native -funroll-loops", "-O3 -march=native -funroll-loops", "-Ofast -march=native -funroll-loops",
    "-O1 -march=native -fdisable-tree-cunrolli", "-O2 -march=native -fdisable-tree-cunrolli", "-O3 -march=native -fdisable-tree-cunrolli", "-Ofast -march=native -fdisable-tree-cunrolli"
]

cpp_files = [
    'eigenrand_real.cpp',
    'pcg_real.cpp',
    'xoshiro_real.cpp',
    'mersenne_twister_real.cpp',
    'std_real.cpp'
]

@rfm.simple_test
class RNGBenchmarkTest(rfm.RegressionTest):
    valid_systems = ['*']
    valid_prog_environs = ['gcc']
    build_system = 'SingleSource'
    sourcesdir = 'src/reframe'
    flags = parameter(flags)
    cpp_file = parameter(cpp_files)

    def __init__(self, flags=None, cpp_file=None):
        super().__init__()

        # Default values for flags and cpp_file if not provided
        self.flags = flags if flags else '-O3 -march=native'
        self.cpp_file = cpp_file if cpp_file else 'std_real.cpp'

        self.valid_systems = ['*']
        self.valid_prog_environs = ['*']
        self.build_system = 'SingleSource'
        self.build_system.cxx = 'g++'
        self.build_system.cxxflags = [self.flags]
        self.sourcepath = self.cpp_file
        # self.sanity_patterns = sn.assert_found(r'PASSED', self.stdout)

        self.tags = {'production'}

    @property
    def name(self):
        return f'{self.__class__.__name__}_flags_{self.flags.replace(" ", "_")}_{os.path.splitext(self.cpp_file)[0]}'

    @run_before('run')
    def set_prerun_cmds(self):
        if self.cpp_file == 'eigenrand_real.cpp':
            build_cmd = f'g++ -Wno-enum-compare {self.flags} -I{os.path.abspath("../extlibs/EigenRand")} -o {os.path.splitext(self.cpp_file)[0]} {self.cpp_file}'
            self.prerun_cmds = [
                f'echo "Building {self.cpp_file} with {self.flags}"',
                'echo "Executing build command: ' + build_cmd + '"',
                build_cmd,
                f'./{os.path.splitext(self.cpp_file)[0]}'
            ]
        if self.cpp_file == 'pcg_real.cpp':
            build_cmd = f'g++ -Wno-enum-compare {self.flags} -I{os.path.abspath("../extlibs/pcg-cpp/include")} -o {os.path.splitext(self.cpp_file)[0]} {self.cpp_file}'
            self.prerun_cmds = [
                f'echo "Building {self.cpp_file} with {self.flags}"',
                'echo "Executing build command: ' + build_cmd + '"',
                build_cmd,
                f'./{os.path.splitext(self.cpp_file)[0]}'
            ]
        if self.cpp_file == 'xoshiro_real.cpp':
            build_cmd = f'g++ -Wno-enum-compare {self.flags} -I{os.path.abspath("../extlibs/Xoshiro-cpp/include")} -o {os.path.splitext(self.cpp_file)[0]} {self.cpp_file}'
            self.prerun_cmds = [
                f'echo "Building {self.cpp_file} with {self.flags}"',
                'echo "Executing build command: ' + build_cmd + '"',
                build_cmd,
                f'./{os.path.splitext(self.cpp_file)[0]}'
            ]
        if self.cpp_file == 'mersenne_twister_real.cpp':
            build_cmd = f'g++ -Wno-enum-compare {self.flags} -o {os.path.splitext(self.cpp_file)[0]} {self.cpp_file}'
            self.prerun_cmds = [
                f'echo "Building {self.cpp_file} with {self.flags}"',
                'echo "Executing build command: ' + build_cmd + '"',
                build_cmd,
                f'./{os.path.splitext(self.cpp_file)[0]}'
            ]
        if self.cpp_file == 'std_real.cpp':
            build_cmd = f'g++ -Wno-enum-compare {self.flags} -o {os.path.splitext(self.cpp_file)[0]} {self.cpp_file}'
            self.prerun_cmds = [
                f'echo "Building {self.cpp_file} with {self.flags}"',
                build_cmd,
                f'./{os.path.splitext(self.cpp_file)[0]}'
            ]

    @performance_function('μs')
    def extract_runtime(self):
        regex_map = {
            'eigenrand_real.cpp': 'EigenRand,(\S+)',
            'pcg_rela.cpp': 'pcg-cpp,(\S+)',
            'xoshiro_real.cpp': 'xoshiro-cpp,(\S+)',
            'mersenne_twister_real.cpp': 'MersenneTwister,(\S+)',
            'std_real.cpp': 'std::uniform_distribution,(\S+)'
        }
        return sn.extractsingle(regex_map[self.cpp_file], self.stdout, 1, float)


    @sanity_function
    def assert_output(self):
        regex_map = {
            'eigenrand_real.cpp': 'EigenRand,(\S+)',
            'pcg_real.cpp': 'pcg-cpp,(\S+)',
            'xoshiro_real.cpp': 'xoshiro-cpp,(\S+)',
            'mersenne_twister_real.cpp': 'MersenneTwister,(\S+)',
            'std_real.cpp': 'std::uniform_distribution,(\S+)'
        }
        return sn.assert_found(regex_map[self.cpp_file], self.stdout)
