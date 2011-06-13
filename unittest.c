#include "minimath.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#define assert_eq(a,b) do {                                     \
  if( fabs((a) - (b)) > 1e-6 )                                  \
  {                                                             \
    printf("Test failed on line %d\n", __LINE__);               \
    return 0;                                                   \
  }                                                             \
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

  double sumref[] = {1.205217, 1.3236178,0.57063374, 1.1164861, 1.0300605};
  double difref[] = {-0.1957338,0.28181938,0.25827212,-0.82205425,-0.076729813};
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

  printf("all tests pass!\n");

  return 0;
}
