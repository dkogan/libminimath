#include "minimath.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#define assert_eq(a,b) do {                                     \
  if( fabs((a) - (b)) > 1e-6 )                                  \
  {                                                             \
    printf("Test failed on line %d\n", __LINE__);               \
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

  // now some matrix inversions
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

  printf("all tests pass!\n");

  return 0;
}
