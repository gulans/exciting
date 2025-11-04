!> @brief Test module for the Fortran implementation of Faddeeva_w.
MODULE faddeeva_test_module
  USE, INTRINSIC :: iso_fortran_env, ONLY: real64
  USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_nan
  USE Faddeeva_module, ONLY: Faddeeva_w

  IMPLICIT NONE
  
  PUBLIC :: run_faddeeva_tests

CONTAINS

  !> @brief Runs a suite of tests for the Faddeeva_w function.
  !> @note Compares results against known values from the original C++ implementation.
  SUBROUTINE run_faddeeva_tests()
    
    INTEGER, PARAMETER :: NTST = 57

    COMPLEX(real64), PARAMETER :: z(NTST) = [ &
        CMPLX(624.2_real64,-0.26123_real64, kind=real64), &
        CMPLX(-0.4_real64,3._real64, kind=real64), &
        CMPLX(0.6_real64,2._real64, kind=real64), &
        CMPLX(-1._real64,1._real64, kind=real64), &
        CMPLX(-1._real64,-9._real64, kind=real64), &
        CMPLX(-1._real64,9._real64, kind=real64), &
        CMPLX(-0.0000000234545_real64,1.1234_real64, kind=real64), &
        CMPLX(-3._real64,5.1_real64, kind=real64), &
        CMPLX(-53._real64,30.1_real64, kind=real64), &
        CMPLX(0.0_real64,0.12345_real64, kind=real64), &
        CMPLX(11._real64,1._real64, kind=real64), &
        CMPLX(-22._real64,-2._real64, kind=real64), &
        CMPLX(9._real64,-28._real64, kind=real64), &
        CMPLX(21._real64,-33._real64, kind=real64), &
        CMPLX(1e5_real64,1e5_real64, kind=real64), &
        CMPLX(1e14_real64,1e14_real64, kind=real64), &
        CMPLX(-3001._real64,-1000._real64, kind=real64), &
        CMPLX(1e160_real64,-1e159_real64, kind=real64), &
        CMPLX(-6.01_real64,0.01_real64, kind=real64), &
        CMPLX(-0.7_real64,-0.7_real64, kind=real64), &
        CMPLX(2.611780000000000e+01_real64, 4.540909610972489e+03_real64, kind=real64), &
        CMPLX(0.8e7_real64,0.3e7_real64, kind=real64), &
        CMPLX(-20._real64,-19.8081_real64, kind=real64), &
        CMPLX(1e-16_real64,-1.1e-16_real64, kind=real64), &
        CMPLX(2.3e-8_real64,1.3e-8_real64, kind=real64), &
        CMPLX(6.3_real64,-1e-13_real64, kind=real64), &
        CMPLX(6.3_real64,1e-20_real64, kind=real64), &
        CMPLX(1e-20_real64,6.3_real64, kind=real64), &
        CMPLX(1e-20_real64,16.3_real64, kind=real64), &
        CMPLX(9._real64,1e-300_real64, kind=real64), &
        CMPLX(6.01_real64,0.11_real64, kind=real64), &
        CMPLX(8.01_real64,1.01e-10_real64, kind=real64), &
        CMPLX(28.01_real64,1e-300_real64, kind=real64), &
        CMPLX(10.01_real64,1e-200_real64, kind=real64), &
        CMPLX(10.01_real64,-1e-200_real64, kind=real64), &
        CMPLX(10.01_real64,0.99e-10_real64, kind=real64), &
        CMPLX(10.01_real64,-0.99e-10_real64, kind=real64), &
        CMPLX(1e-20_real64,7.01_real64, kind=real64), &
        CMPLX(-1._real64,7.01_real64, kind=real64), &
        CMPLX(5.99_real64,7.01_real64, kind=real64), &
        CMPLX(1._real64,0._real64, kind=real64), &
        CMPLX(55._real64,0._real64, kind=real64), &
        CMPLX(-0.1_real64,0._real64, kind=real64), &
        CMPLX(1e-20_real64,0._real64, kind=real64), &
        CMPLX(0._real64,5e-14_real64, kind=real64), &
        CMPLX(0._real64,51._real64, kind=real64), &
        CMPLX(HUGE(0.0_real64),0._real64, kind=real64), &
        CMPLX(-HUGE(0.0_real64),0._real64, kind=real64), &
        CMPLX(0._real64,HUGE(0.0_real64), kind=real64), &
        CMPLX(0._real64,-HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64),HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64),-HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64), & ! NaN equivalent
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), 0.0_real64, kind=real64), &
        CMPLX(0.0_real64, HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64) ]

    COMPLEX(real64), PARAMETER :: w_expected(NTST) = [ &
        CMPLX(-3.7827024551898050745e-7_real64, 0.00090386127643317205733_real64, kind=real64), &
        CMPLX(0.17649062270048168473_real64, -0.021465505394684576168_real64, kind=real64), &
        CMPLX(0.24102507157726921461_real64, 0.060875796634280897459_real64, kind=real64), &
        CMPLX(0.30474420525691259246_real64, -0.20821893820283162729_real64, kind=real64), &
        CMPLX(7.3171310689723780969e34_real64, 8.3218734997144027772e34_real64, kind=real64), &
        CMPLX(0.061569850723632368552_real64, -0.0067600578371657501307_real64, kind=real64), &
        CMPLX(0.39607930076998749190_real64, -5.5931522591166449205e-9_real64, kind=real64), &
        CMPLX(0.082171992267394479433_real64, -0.047012910876436098910_real64, kind=real64), &
        CMPLX(0.0045724600035028164095_real64, -0.0080490079141169182182_real64, kind=real64), &
        CMPLX(0.87463428596080526661_real64, 0.0_real64, kind=real64), &
        CMPLX(0.0046819016496544417437_real64, 0.051073556390130619799_real64, kind=real64), &
        CMPLX(-0.0023193175200187620902_real64, -0.025460054739731556005_real64, kind=real64), &
        CMPLX(9.1146336840563717466e304_real64, 3.9710180714526333377e305_real64, kind=real64), &
        CMPLX(-4.4927207857715598976e281_real64, -2.8019591213423077494e281_real64, kind=real64), &
        CMPLX(2.8209479178093051327e-6_real64, 2.8209479176682577368e-6_real64, kind=real64), &
        CMPLX(2.8209479177387814347e-15_real64, 2.8209479177387814347e-15_real64, kind=real64), &
        CMPLX(-0.000056385128969624435_real64, -0.0001692117551268121746_real64, kind=real64), &
        CMPLX(-5.5860354806708543262e-162_real64, 5.5860354806708543262e-161_real64, kind=real64), &
        CMPLX(0.00016318325137140451888_real64, -0.095232456573009287371_real64, kind=real64), &
        CMPLX(0.69504753678406939989_real64, -1.8916411171103639137_real64, kind=real64), &
        CMPLX(0.00012424182696532796566_real64, 7.1459758263201868885e-7_real64, kind=real64), &
        CMPLX(2.3185873296483533186e-8_real64, 6.1828995457288574857e-8_real64, kind=real64), &
        CMPLX(-0.013342687724350602205_real64, -0.014808709714322076949_real64, kind=real64), &
        CMPLX(1.0_real64, 1.1283791670955127939e-16_real64, kind=real64), &
        CMPLX(0.99999998533107046776_real64, 2.5952720245196788819e-8_real64, kind=real64), &
        CMPLX(-1.4731421795638279504e-15_real64, 0.090727659684127365236_real64, kind=real64), &
        CMPLX(5.7924607788441028458e-18_real64, 0.090727659684127365236_real64, kind=real64), &
        CMPLX(0.088465899352852195347_real64, 1.3708835249574912528e-22_real64, kind=real64), &
        CMPLX(0.034548084541919042437_real64, 2.1116110289517904497e-23_real64, kind=real64), &
        CMPLX(6.6396771995807344007e-36_real64, 0.063082090059258286371_real64, kind=real64), &
        CMPLX(0.0017943523320870264489_real64, 0.095198381480527064794_real64, kind=real64), &
        CMPLX(9.0976037710209799992e-13_real64, 0.070997921072513855099_real64, kind=real64), &
        CMPLX(7.2049510279742166460e-304_real64, 0.020155295647952695387_real64, kind=real64), &
        CMPLX(3.0454360465225073419e-44_real64, 0.056648165176067504293_real64, kind=real64), &
        CMPLX(3.0454360465225073419e-44_real64, 0.056648165176067504293_real64, kind=real64), &
        CMPLX(0.56599287320652734293e-12_real64, 0.056648165176067504293_real64, kind=real64), &
        CMPLX(-0.5659928732065273429e-12_real64, 0.05664816517606750429_real64, kind=real64), &
        CMPLX(0.079688425172165221569_real64, 1.1147446181756167502e-22_real64, kind=real64), &
        CMPLX(0.078171958212473574585_real64, -0.010939136701035766908_real64, kind=real64), &
        CMPLX(0.046700329809904499128_real64, 0.039440389619335341376_real64, kind=real64), &
        CMPLX(0.36787944117144232160_real64, 0.60715770584139372912_real64, kind=real64), &
        CMPLX(0.0_real64, 0.01025968880553683098_real64, kind=real64), &
        CMPLX(0.99004983374916805357_real64, -0.11208866436449538037_real64, kind=real64), &
        CMPLX(0.9999999999999999999999999999999999999999_real64, 1.1283791670955125739e-20_real64, kind=real64), &
        CMPLX(0.99999999999994358104_real64, 0.0_real64, kind=real64), &
        CMPLX(0.011060415485327720154_real64, 0.0_real64, kind=real64), &
        CMPLX(0.0_real64, 0.0_real64, kind=real64), &
        CMPLX(0.0_real64, 0.0_real64, kind=real64), &
        CMPLX(0.0_real64, 0.0_real64, kind=real64), &
        CMPLX(HUGE(0.0_real64), 0.0_real64, kind=real64), &
        CMPLX(0.0_real64, 0.0_real64, kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), 0.0_real64, kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64), &
        CMPLX(HUGE(0.0_real64)/HUGE(0.0_real64), HUGE(0.0_real64)/HUGE(0.0_real64), kind=real64) ]

    COMPLEX(real64) :: fw
    REAL(real64)    :: re_err, im_err, errmax
    INTEGER         :: i

    errmax = 0.0_real64
    WRITE(*,'(A)') "############# w(z) tests #############"

    DO i = 1, NTST
        fw = Faddeeva_w(z(i), 0.0_real64)
        re_err = relerr(REAL(w_expected(i), kind=real64), REAL(fw, kind=real64))
        im_err = relerr(AIMAG(w_expected(i)), AIMAG(fw))
        
        WRITE(*,'(A,ES11.4,"+",ES11.4,"i) = ",ES15.5,"+",ES15.5,"i (vs. ",ES15.5,"+",ES15.5,"i), re/im rel. err. = ",ES11.3,"/",ES11.3)') &
        "w(", REAL(z(i)), AIMAG(z(i)), REAL(fw), AIMAG(fw), &
        REAL(w_expected(i)), AIMAG(w_expected(i)), re_err, im_err
        
        IF (.NOT. ieee_is_nan(re_err)) errmax = MAX(errmax, re_err)
        IF (.NOT. ieee_is_nan(im_err)) errmax = MAX(errmax, im_err)
    END DO

    IF (errmax > 1.0e-13_real64) THEN
        WRITE(*,'(A,ES10.3,A)') "FAILURE -- relative error ", errmax, " too large!"
    ELSE
        WRITE(*,'(A,ES10.3,A)') "SUCCESS (max relative error = ", errmax, ")"
    END IF

  END SUBROUTINE run_faddeeva_tests

  !> @brief Calculates relative error, handling special values.
  FUNCTION relerr(a, b) RESULT(err)
    REAL(real64), INTENT(IN) :: a, b
    REAL(real64) :: err
    IF (ieee_is_nan(a) .OR. ieee_is_nan(b)) THEN
        IF (ieee_is_nan(a) .AND. ieee_is_nan(b)) THEN
            err = 0.0_real64
        ELSE
            err = HUGE(0.0_real64)
        END IF
        RETURN
    END IF
    IF (a == 0.0_real64) THEN
      IF (b == 0.0_real64) THEN
        err = 0.0_real64
      ELSE
        err = HUGE(0.0_real64)
      END IF
    ELSE
      err = ABS((b - a) / a)
    END IF
  END FUNCTION relerr

END MODULE faddeeva_test_module
