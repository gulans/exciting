"""
 Failing and hanging tests.

 Each skipped test should have a TODO(Name) assigned to it here, which
 is also logged on the gitlab issue tracker.

 Each list entry should be a dictionary.

TODO(A/B/H) Issue 51. Review tolerances for all tests migrated from the old test suite (prior to Jan 2021)
TODO(Alex) Convert this into the config file
"""
from src.runner.profile import CompilerBuild, Compiler, Build_type


# Tests that SOMETIMES fail the CI for no discernible reason are considered flakey rather than failing.
#
# Mimics the ctest behaviour: -repeat until-pass:2.
repeat_tests = []

failing_tests = [
    # TODO(Alex) Issue #103. App Tests: GW QP Energies Inconsistent with Reference
    # TODO(Alex) Issue #35. MPI GW calculations do not produce EPS00_GW.OUT
    # I suspect this occurs due to how the GW code distributes MPI processes.
    # If more ranks than q-points, ranks can lie idle starting from rank 0 (opposite from what one usually expects).
    # Also appear to get a failure for Intel 2019 serial (only occurs on the CI)
    {'name': 'GW/LDA_PW-gw-Si',
     'comment': 'Serial calcs give erroneous QP energies w.r.t. reference.'
                'MPI GW calculations do not produce EPS00_GW.OUT',
     'tags': [CompilerBuild(Compiler.all, Build_type.serial),
              CompilerBuild(Compiler.all, Build_type.mpiandsmp)]
     },

    # TODO(Sebastian) Issue 40
    {'name': 'properties/LDA_PZ-wannier-SiC',
     'comment': '',
     'tags': [CompilerBuild(Compiler.all, Build_type.all)]
     },

    # TODO(Cecilia) Issue 54. HSE and PBE0 sometimes fail in our CI pipeline within Intel MPIANDSMP build
    # Repeating them changes the results each time but the diff is never < tolerances.
    # Only occurs in the CI pipeline.
    {'name': 'hybrid/HSE-Si',
     'comment': 'Test is flakey when run in the CI with Intel parallel build.',
     'tags': [CompilerBuild(Compiler.intel, Build_type.mpiandsmp)]
     },

    # TODO(Cecilia) Issue 54
    {'name': 'hybrid/PBE0-Si',
     'comment': 'Test is flakey when run in the CI with Intel parallel build.',
     'tags': [CompilerBuild(Compiler.intel, Build_type.mpiandsmp)]
     },

    # TODO(Alex) Issue 111.
    {'name': 'properties/LDA_PW-optical_properties-Si',
     'comment': 'Epsilon 11 and 33 do not agree with serial reference values. \n'
                '  - CHI_111.OUT also reviews reviewing.',
     'tags': [CompilerBuild(Compiler.all, Build_type.all)]
     },

    # TODO(Alex) Issue 112
    {'name': 'properties/LDA_PW-ldos-Si',
     'comment': 'Numbers are completely different - reference data needs reviewing.',
     'tags': [CompilerBuild(Compiler.all, Build_type.all)]
     },

    # TODO(Maria) Issue 55
    {'name': 'properties/LDA_PW-transport-Si',
     'comment': 'Test is flakey when run in the CI with GCC builds: Test outputs are not written.\n'
                '  Results with Intel for THERMALCOND_11 and ELECTCOND_11 have large absolute errors.',
     'tags': [CompilerBuild(Compiler.all, Build_type.all)]
     },

    # TODO(Bene) Issue 113
    {'name': 'properties/LDA_PW-SO-spintexture-Si',
     'comment': 'GCC calculations disagree with reference.',
     'tags': [CompilerBuild(Compiler.all, Build_type.all)]
     },

    # TODO(Sven) Issue #39
    # Andris comment: I have tried several runs with a different number of threads and MPI ranks,
    # but the calculations always converged in 19-20 iterations and always to the same energy with the threshold.
    # Alex: Update. Test gives different accuracy each time it runs in the CI
    # For example, total energy difference 1.34e-06 Ha, seems too high. See MR !173 for details.
    {'name': 'groundstate/LDA_PW-noncollinear-Fe',
     'comment': 'Variation in total energy w.r.t. reference is too large.',
     'tags': [CompilerBuild(Compiler.all, Build_type.all)]
     },

    # TODO(Alex) Issue 101. Also see Issue #75
    {'name': 'groundstate/LDA_PW-collinear-Fe',
     'comment': 'Most energies differ to reference by ~1.e-7. '
                'scl%Sum of eigenvalues by ~ 1.e-6. '
                'DOS at Fermi differs by 5.7e-04. ',
     'tags': [CompilerBuild(Compiler.gcc, Build_type.all),
              CompilerBuild(Compiler.intel, Build_type.all)]
     }
]

# These always have to be skipped else the suite hangs 
hanging_tests = [
    # TODO(Alex) Issue #36 chargedensityplot hangs when running with np > 1 cores
    # Should remove chargedensityplot from this and add it as a property
    {'name': 'groundstate/PBE-Al',
     'comment': 'chargedensityplot hangs at an allgatherv call.',
     'tags': [CompilerBuild(Compiler.all, Build_type.mpiandsmp)]}
]