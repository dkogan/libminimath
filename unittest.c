#include "minimath.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#define assert_eq(a,b) do {                                     \
  if( fabs((a) - (b)) > 1e-2 )                                  \
  {                                                             \
    printf("Test failed on line %d. Err: %f\n", __LINE__, fabs((a) - (b)));   \
    return 1;                                                   \
  }                                                             \
} while(0)

#define assert_vector_elemeq6(a,b) do {                         \
  for(int i=0; i<6; i++) assert_eq(a,b);                        \
} while(0)

#define assert_vector_elemeq5(a,b) do {                         \
  for(int i=0; i<5; i++) assert_eq(a,b);                        \
} while(0)

#define assert_vector_elemeq3(a,b) do {                         \
  for(int i=0; i<3; i++) assert_eq(a,b);                        \
} while(0)

int main(void)
{
  double a5[] = {0.5047416,0.80271857,0.41445293,0.14721594,0.47666535};
  double b5[] = {0.7004754,0.52089919,0.15618081,0.96927019,0.55339517};

  double v5[5];

  // first, make sure norms and dot-products work
  {
    assert_eq( norm2_vec(5, a5),     1.31977481447686 );
    assert_eq( dot_vec  (5, a5, b5), 1.24290044685416 );
  }

  // now test various flavors of vector addition/subtraction
  {
    add_vec_vout(5, a5, b5, v5);
    assert_vector_elemeq5(a5[i]+b5[i], v5[i]);

    sub_vec_vout(5, a5, b5, v5);
    assert_vector_elemeq5(a5[i]-b5[i], v5[i]);

    memcpy(v5, a5, sizeof(a5));
    add_vec(5, v5, b5);
    assert_vector_elemeq5(a5[i]+b5[i], v5[i]);

    memcpy(v5, a5, sizeof(a5));
    sub_vec(5, v5, b5);
    assert_vector_elemeq5(a5[i]-b5[i], v5[i]);

    memcpy(v5, a5, sizeof(a5));
    add_vec_vaccum(5, a5, b5, v5);
    assert_vector_elemeq5(2.0*a5[i]+b5[i], v5[i]);

    memcpy(v5, a5, sizeof(a5));
    sub_vec_vaccum(5, a5, b5, v5);
    assert_vector_elemeq5(2.0*a5[i]-b5[i], v5[i]);

    // now the scaled versions
    add_vec_vout_scaled(5, a5, b5, v5, 3.0);
    assert_vector_elemeq5(a5[i]+3.0*b5[i], v5[i]);

    memcpy(v5, a5, sizeof(a5));
    sub_vec_scaled(5, v5, b5, 3.0);
    assert_vector_elemeq5(a5[i]-3.0*b5[i], v5[i]);
  }

  // now do various matrix things, so define some test matrices
  // symmetric matrix:
  // [0.41974018,0.92003082,  1.004335,   1.30167,0.87259315],
  // [0.92003082,0.91868919, 1.1176797, 1.2953584,0.60691799],
  // [  1.004335, 1.1176797, 1.3854073, 1.0001464,0.73574336],
  // [   1.30167, 1.2953584, 1.0001464, 1.8231297, 1.1258482],
  // [0.87259315,0.60691799,0.73574336, 1.1258482,0.40574078]
  double s5[] = {0.41974018,0.92003082,  1.004335,   1.30167,0.87259315,
                            0.91868919, 1.1176797, 1.2953584,0.60691799,
                                        1.3854073, 1.0001464,0.73574336,
                                                   1.8231297, 1.1258482, 
                                                              0.40574078};
  double a_x_s[] = {1.9741972, 2.1450465, 2.4762389, 2.9163754, 1.5916948};

  double m35[] =
    { 0.27385016, 0.37604514, 0.91209475, 0.93284312, 0.27634535,
      0.79741259, 0.94574321, 0.23214214, 0.27074102, 0.44384911,
      0.68883162,0.069744802,  0.2323944,0.073613799, 0.81272502};

  double m53[] =
    {
      0.42868455, 0.96784151, 0.50407678,
      0.50540621, 0.22897339, 0.15036613,
      0.93139262, 0.22178075, 0.78665219,
      0.57617856, 0.87582742, 0.20284634,
      0.60838134,0.037992089,0.060601891};

  double a_x_m53 [] = {1.3829093,0.91127402,0.75990955};
  double a_x_m35t[] = {1.087156, 1.5092898,0.89821894};
  double v3[3];
  double v6[6];

  // for 3-way symmetric multiplication
  // pdl> p $a
  // [
  //  [ 1.7154671,0.32440803, 1.5281059],
  //  [0.32440803, 1.3319233, 1.4894104],
  //  [ 1.5281059, 1.4894104, 1.8696738]
  // ]
  double s3_a[] = {1.7154671,0.32440803, 1.5281059, 1.3319233, 1.4894104, 1.8696738};

  // pdl> p $b
  // [
  //  [ 1.3211548, 1.6953079,0.46298047],
  //  [ 1.6953079, 1.4032645, 1.2392179],
  //  [0.46298047, 1.2392179, 1.5834584]
  // ]
  double s3_b[] = {1.3211548, 1.6953079,0.46298047, 1.4032645, 1.2392179, 1.5834584};

  // pdl> p $a x $b x $a
  // [
  //  [ 13.276035, 13.530868, 19.975459],
  //  [ 13.530868, 12.970227, 19.287342],
  //  [ 19.975459, 19.287342, 28.997443]
  // ]
  double s3_aba[] = {13.276035, 13.530868, 19.975459, 12.970227, 19.287342, 28.997443};

  // pdl> p $a3 = random(3)
  // [0.93668206, 0.5618703,0.71166218]
  // pdl> p $b3 = random(3)
  // [0.86716519,0.22144078,0.83853122]
  // pdl> p ($a3 x $a x $b3->transpose)
  // [
  //  [ 5.9799184]
  // ]
  double a3[]        = {0.93668206, 0.5618703,0.71166218};
  double b3[]        = {0.86716519,0.22144078,0.83853122};
  double conj_result = 5.9799184;

  // symmetric multiplication
  {
    mul_vec5_sym55_vout(a5, s5, v5);
    assert_vector_elemeq5(v5[i], a_x_s[i]);
    mul_vec5_sym55_vout_scaled(a5, s5, v5, -3.0);
    assert_vector_elemeq5(v5[i], -3.0*a_x_s[i]);

    memcpy(v5, a5, sizeof(a5));
    mul_vec5_sym55(v5, s5);
    assert_vector_elemeq5(v5[i], a_x_s[i]);
    memcpy(v5, a5, sizeof(a5));
    mul_vec5_sym55_scaled(v5, s5, -3.0);
    assert_vector_elemeq5(v5[i], -3.0*a_x_s[i]);

    memcpy(v5, a5, sizeof(a5));
    mul_vec5_sym55_vaccum(a5, s5, v5);
    assert_vector_elemeq5(v5[i], a_x_s[i] + a5[i]);
    memcpy(v5, a5, sizeof(a5));
    mul_vec5_sym55_vaccum_scaled(a5, s5, v5, -3.0);
    assert_vector_elemeq5(v5[i], -3.0*a_x_s[i] + a5[i]);

    // 3-way symmetric multiplication
    mul_sym33_sym33_sym33_vout(s3_a, s3_b, v6);
    assert_vector_elemeq6(v6[i], s3_aba[i]);

    // conjugation
    assert_eq( conj_3(a3, s3_a, b3),
               conj_result );
  }

  // general multiplication
  {
    mul_vec5_gen53_vout(a5, m53, v3);
    assert_vector_elemeq3(v3[i], a_x_m53[i]);
    mul_vec5_gen53_vout_scaled(a5, m53, v3, -3.0);
    assert_vector_elemeq3(v3[i], -3.0*a_x_m53[i]);

    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen53(v5, m53);
    assert_vector_elemeq3(v5[i], a_x_m53[i]);
    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen53_scaled(v5, m53, -3.0);
    assert_vector_elemeq3(v5[i], -3.0*a_x_m53[i]);

    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen53_vaccum(a5, m53, v5);
    assert_vector_elemeq3(v5[i], a_x_m53[i] + a5[i]);
    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen53_vaccum_scaled(a5, m53, v5, -3.0);
    assert_vector_elemeq3(v5[i], -3.0*a_x_m53[i] + a5[i]);
  }

  // general transposed multiplication
  {
    mul_vec5_gen35t_vout(a5, m35, v3);
    assert_vector_elemeq3(v3[i], a_x_m35t[i]);
    mul_vec5_gen35t_vout_scaled(a5, m35, v3, -3.0);
    assert_vector_elemeq3(v3[i], -3.0*a_x_m35t[i]);

    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen35t(v5, m35);
    assert_vector_elemeq3(v5[i], a_x_m35t[i]);
    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen35t_scaled(v5, m35, -3.0);
    assert_vector_elemeq3(v5[i], -3.0*a_x_m35t[i]);

    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen35t_vaccum(a5, m35, v5);
    assert_vector_elemeq3(v5[i], a_x_m35t[i] + a5[i]);
    memcpy(v5, a5, sizeof(a5));
    mul_vec5_gen35t_vaccum_scaled(a5, m35, v5, -3.0);
    assert_vector_elemeq3(v5[i], -3.0*a_x_m35t[i] + a5[i]);
  }

  // now some matrix inversions. Symmetric
  {
    double s3[] =
      { 0.471011, 1.6661985 , 0.98615889,
                  0.32707543, 1.0342404,
                              0.49936779};

    double cofactors[6];
    double m[9];
    double det = cofactors_sym3(s3, cofactors);

    assert_eq(det, 1.26747089766342);

    // multiplied together, these should be identity*determinant
    mul_sym33_sym33_scaled_out(s3, cofactors, m, 1.0/det);

    for(int i=0; i<9; i++)
    {
      assert_eq(m[i], i%4 ? 0 : 1);
    }

  }

  // upper/lower triangular inversions.
  {
      /* Comes from this python session:

In [13]: n = 2; m = np.triu(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.57218581  0.96981434]
 [ 0.          0.44706582]]
[[ 1.74768404 -3.7912293 ]
 [ 0.          2.23680709]]

In [14]: n = 3; m = np.triu(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.3143549   0.48010021  0.98662155]
 [ 0.          0.17672097  0.241954  ]
 [ 0.          0.          0.71325572]]
[[ 3.18111791 -8.64218538 -1.46868528]
 [ 0.          5.65863807 -1.91955014]
 [ 0.          0.          1.40202171]]

In [15]: n = 4; m = np.triu(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.69412981  0.77539124  0.55602195  0.97512905]
 [ 0.          0.38083716  0.05308887  0.53878659]
 [ 0.          0.          0.65848863  0.11064253]
 [ 0.          0.          0.          0.85354224]]
[[ 1.44065272 -2.93319458 -0.97999346  0.33269888]
 [ 0.          2.62579416 -0.21169756 -1.63005399]
 [ 0.          0.          1.5186291  -0.19685606]
 [ 0.          0.          0.          1.17158817]]

In [16]: n = 5; m = np.triu(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.86692312  0.0431848   0.87871003  0.24380479  0.5625632 ]
 [ 0.          0.31158425  0.32071057  0.98432269  0.6176679 ]
 [ 0.          0.          0.29743654  0.03104779  0.94812118]
 [ 0.          0.          0.          0.89684428  0.53583566]
 [ 0.          0.          0.          0.          0.55817574]]
[[ 1.15350482 -0.1598729  -3.23539034 -0.02610461  4.53505711]
 [ 0.          3.20940484 -3.46053666 -3.40265092  5.59308334]
 [ 0.          0.          3.36206165 -0.11639096 -5.59908856]
 [ 0.          0.          0.          1.11502077 -1.07039387]
 [ 0.          0.          0.          0.          1.79155045]]

In [17]: n = 2; m = np.tril(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.62669575  0.        ]
 [ 0.16808425  0.05153309]]
[[  1.59567063   0.        ]
 [ -5.20456104  19.40500797]]

In [18]: n = 3; m = np.tril(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.73281257  0.          0.        ]
 [ 0.15663661  0.6252232   0.        ]
 [ 0.32331516  0.49048184  0.68948659]]
[[ 1.36460541  0.          0.        ]
 [-0.34187337  1.59942882  0.        ]
 [-0.39669363 -1.13778978  1.45035454]]

In [19]: n = 4; m = np.tril(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.16449536  0.          0.          0.        ]
 [ 0.38291535  0.79820702  0.          0.        ]
 [ 0.2376585   0.53143921  0.41291111  0.        ]
 [ 0.55820531  0.98338752  0.4130384   0.59041151]]
[[  6.07919888e+00   3.97783041e-16  -3.37463328e-16  -4.41627891e-32]
 [ -2.91630931e+00   1.25280783e+00   1.91556043e-16  -1.39089610e-16]
 [  2.54456211e-01  -1.61243225e+00   2.42182873e+00   3.31153192e-16]
 [ -1.06820260e+00  -9.58651938e-01  -1.69425606e+00   1.69373393e+00]]

In [20]: n = 5; m = np.tril(np.random.rand(n,n)); print m; print np.linalg.inv(m)
[[ 0.00096286  0.          0.          0.          0.        ]
 [ 0.02479833  0.76528692  0.          0.          0.        ]
 [ 0.58960537  0.90182694  0.15462504  0.          0.        ]
 [ 0.64151936  0.42728706  0.50315051  0.86653404  0.        ]
 [ 0.76332981  0.1009956   0.25110659  0.21694558  0.86720055]]
[[  1.03857381e+03  -2.54528288e-16   2.90889472e-16  -7.71800583e-17
    9.98482526e-16]
 [ -3.36539136e+01   1.30669945e+00   0.00000000e+00   3.36914461e-17
   -3.00309860e-17]
 [ -3.76393566e+03  -7.62112526e+00   6.46725801e+00   0.00000000e+00
   -3.79560274e-15]
 [  1.43322799e+03   3.78085244e+00  -3.75519490e+00   1.15402276e+00
    1.40581191e-15]
 [ -1.78919266e+02   1.10874542e+00  -9.33230681e-01  -2.88699237e-01
    1.15313580e+00]]

      */
      {
          // upper
          {
              double s[]       = {0.57218581,0.96981434,0.44706582};
              double inv_ref[] = {1.74768404,-3.7912293,2.23680709};
              double cofactors[sizeof(s)/sizeof(s[0])];
              double det       = cofactors_ut2(s, cofactors);
              // I make sure that inverse = cofactors/det
              for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                  assert_eq(inv_ref[i], cofactors[i]/det);
          }
          {
              double s[]       = {0.3143549,0.48010021,0.98662155,0.17672097,0.241954,0.71325572};
              double inv_ref[] = {3.18111791,-8.64218538,-1.46868528,5.65863807,-1.91955014,1.40202171};
              double cofactors[sizeof(s)/sizeof(s[0])];
              double det       = cofactors_ut3(s, cofactors);
              // I make sure that inverse = cofactors/det
              for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                  assert_eq(inv_ref[i], cofactors[i]/det);
          }
          {
              double s[]       = {0.69412981,0.77539124,0.55602195,0.97512905,0.38083716,0.05308887,0.53878659,0.65848863,0.11064253,0.85354224};
              double inv_ref[] = {1.44065272,-2.93319458,-0.97999346,0.33269888,2.62579416,-0.21169756,-1.63005399,1.5186291,-0.19685606,1.17158817};
              double cofactors[sizeof(s)/sizeof(s[0])];
              double det       = cofactors_ut4(s, cofactors);
              // I make sure that inverse = cofactors/det
              for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                  assert_eq(inv_ref[i], cofactors[i]/det);
          }
          {
              double s[]       = {0.86692312,0.0431848,0.87871003,0.24380479,0.5625632,0.31158425,0.32071057,0.98432269,0.6176679,0.29743654,0.03104779,0.94812118,0.89684428,0.53583566,0.55817574};
              double inv_ref[] = {1.15350482,-0.1598729,-3.23539034,-0.02610461,4.53505711,3.20940484,-3.46053666,-3.40265092,5.59308334,3.36206165,-0.11639096,-5.59908856,1.11502077,-1.07039387,1.79155045};
              double cofactors[sizeof(s)/sizeof(s[0])];
              double det       = cofactors_ut5(s, cofactors);
              // I make sure that inverse = cofactors/det
              for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                  assert_eq(inv_ref[i], cofactors[i]/det);
          }

          // lower
          {
              {
                  double s[]       = {0.62669575,0.16808425,0.05153309};
                  double inv_ref[] = {1.59567063,-5.20456104,19.40500797};
                  double cofactors[sizeof(s)/sizeof(s[0])];
                  double det       = cofactors_lt2(s, cofactors);
                  // I make sure that inverse = cofactors/det
                  for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                      assert_eq(inv_ref[i], cofactors[i]/det);
              }
              {
                  double s[]       = {0.73281257,0.15663661,0.6252232,0.32331516,0.49048184,0.68948659};
                  double inv_ref[] = {1.36460541,-0.34187337,1.59942882,-0.39669363,-1.13778978,1.45035454};
                  double cofactors[sizeof(s)/sizeof(s[0])];
                  double det       = cofactors_lt3(s, cofactors);
                  // I make sure that inverse = cofactors/det
                  for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                      assert_eq(inv_ref[i], cofactors[i]/det);
              }
              {
                  double s[]       = {0.16449536,0.38291535,0.79820702,0.2376585,0.53143921,0.41291111,0.55820531,0.98338752,0.4130384,0.59041151};
                  double inv_ref[] = {6.07919888e+00,-2.91630931e+00,1.25280783e+00,2.54456211e-01,-1.61243225e+00,2.42182873e+00,-1.06820260e+00,-9.58651938e-01,-1.69425606e+00,1.69373393e+00};
                  double cofactors[sizeof(s)/sizeof(s[0])];
                  double det       = cofactors_lt4(s, cofactors);
                  // I make sure that inverse = cofactors/det
                  for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                      assert_eq(inv_ref[i], cofactors[i]/det);
              }

              {
                  double s[]       = {0.00096286,0.02479833,0.76528692,0.58960537,0.90182694,0.15462504,0.64151936,0.42728706,0.50315051,0.86653404,0.76332981,0.1009956,0.25110659,0.21694558,0.86720055};
                  double inv_ref[] = {1.03857381e+03,-3.36539136e+01,1.30669945e+00,-3.76393566e+03,-7.62112526e+00,6.46725801e+00,1.43322799e+03,3.78085244e+00,-3.75519490e+00,1.15402276e+00,-1.78919266e+02,1.10874542e+00,-9.33230681e-01,-2.88699237e-01,1.15313580e+00};
                  double cofactors[sizeof(s)/sizeof(s[0])];
                  double det       = cofactors_lt5(s, cofactors);
                  // I make sure that inverse = cofactors/det
                  for(int i=0; i<sizeof(s)/sizeof(s[0]); i++)
                      assert_eq(inv_ref[i], cofactors[i]/det);
              }
          }
      }
  }
  printf("all tests pass!\n");

  return 0;
}

