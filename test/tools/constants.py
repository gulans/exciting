"""
 Environment and directory variables
"""

from collections import namedtuple

from tools.utils import Build_type, build_type_enum_to_str

_binary_names    = [binary for binary in build_type_enum_to_str.values()]
_exe_ref         = Build_type.puresmp
_default_np      = {Build_type.serial: 1, Build_type.puresmp: 1, Build_type.purempi: 2, Build_type.mpiandsmp: 2}
_default_threads = {Build_type.serial: 1, Build_type.puresmp: 2, Build_type.purempi: 1, Build_type.mpiandsmp: 2}
_action_choices  = ['run', 'ref', 'clean']

Defaults = namedtuple('Defaults', ['max_time',       # Time after which a test is killed (in seconds)
                                   'test_farm',      # Test farm directory
                                   'species',        # Species file directory
                                   'input_file',     # Input file for exciting
                                   'main_output',    # Main output from exciting
                                   'run_dir',        # Run directory for tests
                                   'ref_dir',        # Directory for test reference data
                                   'exe_dir',        # Location of exciting executable
                                   'init_default',   # Template for init xml
                                   'ignored_output', # Output files to not reference
                                   'action_choices', # Action choices for test script 
                                   'binary_names',   # List of exciting executables
                                   'binary_mpismp',  # Exciting execuatable with smp and mpi parallelisation
                                   'binary_purempi', # Exciting execuatable with mpi parallelisation
                                   'binary_smp',     # Exciting execuatable with smp parallelisation
                                   'binary_serial',  # Serial exciting execuatable
                                   'default_np',     # Dict of default MPI processes per executable
                                   'default_threads',# Dict of default threads per executable 
                                   'exe_ref'         # Executable for running the reference calculations
                                   ])

# Define an immutable instance of the default settings
# Access like settings.max_time
settings = Defaults(max_time       = 1800,
                    test_farm      = 'test_farm',  
                    species        = '../species',
                    input_file     = 'input.xml',
                    main_output    = 'INFO.OUT',    
                    run_dir        = 'run',       
                    ref_dir        = 'ref',        
                    exe_dir        = '../../../../bin/',        
                    init_default   = 'xml/init_templates/init_default.xml' ,   
                    ignored_output = ['STATE.OUT', 'OCC', 'EVEC', 'EVALSV', 'EVALFV', 'APWCMT', 'SYM', 
                                      'PMAT', 'FERMISURF', 'RMSDVEFF', 'LOCMT', 'EXCLI', 'SCCLI'], 
                    action_choices  = _action_choices,
                    binary_names    = _binary_names,
                    binary_mpismp   = Build_type.mpiandsmp,
                    binary_purempi  = Build_type.purempi,
                    binary_smp      = Build_type.puresmp,
                    binary_serial   = Build_type.serial,
                    default_np      = _default_np,
                    default_threads = _default_threads,
                    exe_ref         = _exe_ref
                    )
