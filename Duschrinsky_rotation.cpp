//
// Created by phyzch on 5/8/21.
//
#include "util.h"
#include "system.h"

void detector::compute_Duschrinsky_rotation_overlap(){
    // mode 1 and mode 2 is vibrational coordinate we perform Duschrinsky rotation
    int n_prime, m_prime, n , m;

    int Nonzero_element_upper_cutoff ; // |n' , m'> they have their limit for their maximum nonzero quanta
    int index_x, index_y;

    double omega_y = mfreq[0][2];
    double omega_x = mfreq[0][1];

    double r_plus = double(1)/2 * (sqrt(omega_y / omega_x)  + sqrt(omega_x / omega_y) );
    double r_minus = double(1)/2 * (sqrt(omega_y / omega_x) - sqrt(omega_x / omega_y) );

    // Duschrinsky_rotation_overlap : [Maximum_mode , Maximum_mode, Maximum_mode, Maximum_mode]
    // its element [n',m', n, m] stands for <n,m| n',m'> here n' , m' is rotated Harmonic oscillator basis
    // n' is for coordinate x', m' is for coordinate y'. n is for coordinate x, m is coordinate y.
    int Maximum_nmax = max( nmax[0][1], nmax[0][2] );


    Duschrinsky_rotation_Overlap.resize(Maximum_nmax + 1);
    for(n_prime=0;n_prime <= Maximum_nmax;n_prime++){
        Duschrinsky_rotation_Overlap[ n_prime ].resize(Maximum_nmax + 1);
        for(m_prime=0; m_prime <= Maximum_nmax; m_prime++ ){
            Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ].resize(2 * Maximum_nmax + 1  );
            for(n=0;n <= 2 * Maximum_nmax;n++){
                Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ][ n ].resize(2 * Maximum_nmax + 1);
                for(m =0; m <= 2 * Maximum_nmax; m++ ){
                    Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ][ n ][ m ] = 0 ;
                }
            }
        }
    }

    // same ground state
    Duschrinsky_rotation_Overlap[0][0][0][0] = 1;
    // compute [0 , m', n, m]
    for(m_prime=1; m_prime<= Maximum_nmax; m_prime++){
        // compute Duschrinsky_rotation_Overlap[0][m']:
        Nonzero_element_upper_cutoff =  m_prime - 1;
        // n, m is index for Duschrinsky_rotation_Overlap[0][m'-1]
        // b_{y}'^{+} operator
        for( n = 0; n <= Nonzero_element_upper_cutoff ; n ++){
            for(m = 0; m <= Nonzero_element_upper_cutoff ; m ++){
                // we only have to compute contribution from nonzero coefficient.
                if(Duschrinsky_rotation_Overlap[0][m_prime - 1][n][m] == 0){
                    continue ;
                }

                // cos(\theta) b_{y}^{+} operate upon |0, m'-1> .
                index_x = n ;
                index_y = m + 1;
                if(  (Maximum_nmax >= index_x) and (index_x >= 0) and (Maximum_nmax >= index_y) and (index_y >= 0 )   ){
                    Duschrinsky_rotation_Overlap[0][m_prime][index_x][index_y] =
                            Duschrinsky_rotation_Overlap[0][m_prime][index_x][index_y] +
                            sqrt( double(index_y) ) / sqrt( double(m_prime) ) * cos(rotation_angle) *
                            Duschrinsky_rotation_Overlap[0][m_prime-1][n][m];
                }

                // bx^{+} operator
                index_x = n + 1;
                index_y = m;
                if( (Maximum_nmax >= index_x) and (index_x >= 0) and (Maximum_nmax >= index_y) and (index_y >= 0 )   ){
                    Duschrinsky_rotation_Overlap[0][m_prime][index_x][index_y] =
                            Duschrinsky_rotation_Overlap[0][m_prime][index_x][index_y] +
                            sqrt( double(index_x) ) / sqrt( double(m_prime) ) * (- r_plus * sin(rotation_angle) ) *
                            Duschrinsky_rotation_Overlap[0][m_prime - 1][n][m];
                }

                // bx operator
                index_x = n - 1;
                index_y = m;
                if( (Maximum_nmax >= index_x) and (index_x >= 0) and (Maximum_nmax >= index_y) and (index_y >= 0 )   ){
                    Duschrinsky_rotation_Overlap[0][m_prime][index_x][index_y] =
                            Duschrinsky_rotation_Overlap[0][m_prime][index_x][index_y] +
                            sqrt( double(n) ) / sqrt( double(m_prime) ) * (- r_minus * sin(rotation_angle)) *
                            Duschrinsky_rotation_Overlap[0][m_prime - 1][n][m];
                }

            }
        }

    }

    // compute [n', m', n, m]  n' = 1, 2, 3, etc. m'= 0, 1, 2, 3, etc.
    for(n_prime = 1; n_prime <= Maximum_nmax; n_prime ++ ){
        for(m_prime = 0; m_prime <= Maximum_nmax; m_prime ++ ){
            Nonzero_element_upper_cutoff = n_prime + m_prime -1;
            // n, m is index for Duschrinsky_rotation_Overlap[n' -1][m']
            // b_{x}^{+}' operator
            for(n = 0 ; n <= Nonzero_element_upper_cutoff; n++ ){
                for(m =0; m <= Nonzero_element_upper_cutoff; m++ ){
                    // we only have to compute contribution from nonzero coefficient.
                    if(Duschrinsky_rotation_Overlap[n_prime - 1][m_prime][n][m] == 0){
                        continue ;
                    }

                    // b_{x}^{+} operator
                    index_x = n + 1;
                    index_y = m;
                    if(  (index_x >= 0) and (index_y >= 0 )   ){
                        Duschrinsky_rotation_Overlap[n_prime][m_prime][index_x][index_y] =
                                Duschrinsky_rotation_Overlap[n_prime][m_prime][index_x][index_y] +
                                sqrt( double(index_x) ) / sqrt(  double(n_prime)  ) * cos(rotation_angle) *
                                Duschrinsky_rotation_Overlap[n_prime-1][m_prime][n][m];
                    }

                    // b_{y}^{+} operator
                    index_x = n;
                    index_y = m + 1;
                    if(   (index_x >= 0)  and (index_y >= 0 )   ){
                        Duschrinsky_rotation_Overlap[n_prime][m_prime][index_x][index_y] =
                                Duschrinsky_rotation_Overlap[n_prime][m_prime][index_x][index_y] +
                                sqrt(  double(index_y)  ) / sqrt( double(n_prime) ) * r_plus * sin(rotation_angle) *
                                Duschrinsky_rotation_Overlap[n_prime - 1][m_prime][n][m];
                    }

                    // b_{y} operator
                    index_x = n;
                    index_y = m - 1;
                    if(  (index_x >= 0)  and (index_y >= 0 )   ){
                        Duschrinsky_rotation_Overlap[n_prime][m_prime][index_x][index_y] =
                                Duschrinsky_rotation_Overlap[n_prime][m_prime][index_x] [index_y] +
                                sqrt( double(m) ) / sqrt(double(n_prime)) * (- r_minus) * sin(rotation_angle) *
                                Duschrinsky_rotation_Overlap[n_prime - 1] [m_prime][n][m];
                    }

                }
            }

        }

    }


}

void detector:: compute_shifted_Duschrinsky_rotation_overlap( ){
    // Notice this function 's computational cost grow as N^{8} where N is maximum mode quanta set in x,y coordinate.
    int i , j;

    double omega_y = mfreq[0][2];
    double omega_x = mfreq[0][1];

    double r_plus = double(1)/2 * (sqrt(omega_y / omega_x)  + sqrt(omega_x / omega_y) );
    double r_minus = double(1)/2 * (sqrt(omega_y / omega_x) - sqrt(omega_x / omega_y) );

    double alpha_down = coupling_strength_to_mode0[0] / mfreq[0][1];
    double alpha_up = - coupling_strength_to_mode0[1] / mfreq[0][1];

    double Coeff = 1 / ( pow(cos(rotation_angle) , 2) + r_plus * pow(sin(rotation_angle), 2) );
    double A = Coeff * (alpha_down * cos(rotation_angle) - alpha_up * r_minus * r_plus * pow(sin(rotation_angle), 2 ) ) ;
    double B = Coeff * (- alpha_down * r_plus * sin(rotation_angle) - alpha_up * r_minus * sin(rotation_angle) * cos(rotation_angle)) ;
    double C = alpha_up * cos(rotation_angle) -r_minus * sin(rotation_angle) * Coeff * (alpha_down * r_plus * sin(rotation_angle) + alpha_up * r_minus * sin(rotation_angle) * cos(rotation_angle) );
    double D = alpha_up * r_plus * sin(rotation_angle) - r_minus * sin(rotation_angle) * Coeff * (- alpha_down * cos(rotation_angle) + alpha_up * r_minus * r_plus * pow(sin(rotation_angle) , 2)  ) ;

    double A_prime = A - alpha_up;
    double B_prime = B;
    double C_prime = C - alpha_down;
    double D_prime = D;
    int n_prime, m_prime,  n, m;
    int l, p ,l_prime, p_prime;

    double E  = -cos(rotation_angle) * (A * C + B * D) + r_plus * sin(rotation_angle) * (B * C - A * D);

    double Overlap_coeff = exp ( (- ( pow(alpha_down,2) + pow(alpha_up,2) ) +
                                  alpha_down * alpha_up * cos(rotation_angle) - E  )  / 2 );

    int Maximum_nmax = max( nmax[0][1], nmax[0][2] );
    // This is overlap for shifted Harmonic oscillator <n',m' | n, m>
    shifted_Duschrinsky_rotation_Overlap.resize(Maximum_nmax + 1);
    for(n_prime=0; n_prime<=Maximum_nmax; n_prime++){
        shifted_Duschrinsky_rotation_Overlap[ n_prime ].resize(Maximum_nmax + 1);
        for(m_prime=0; m_prime<= Maximum_nmax; m_prime++ ){
            shifted_Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ].resize(Maximum_nmax + 1);
            for(n=0; n<=Maximum_nmax; n++){
                shifted_Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ][ n ].resize(Maximum_nmax + 1);
                for(m =0; m <= Maximum_nmax; m++ ){
                    shifted_Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ][ n ][ m ] = 0 ;
                }
            }
        }
    }

    // Coeff1 for |n', m'> . Coeff2 for |n,m>
    double Coeff1_x ;
    double Coeff1_y;
    double Coeff2_x ;
    double Coeff2_y;
    double Coeff_prod;

    double n_prime_factorial;
    double n_factorial ;
    double m_prime_factorial ;
    double m_factorial ;

    double l_prime_factorial ;
    double l_factorial ;
    double p_prime_factorial ;
    double p_factorial ;

    double n_prime_minus_l_prime_factorial ;
    double m_prime_minus_p_prime_factorial ;
    double n_minus_l_factorial ;
    double m_minus_p_factorial ;


    for(n_prime = 0; n_prime <= Maximum_nmax; n_prime ++ ){
        n_prime_factorial = factorial(n_prime);
        for(m_prime = 0; m_prime <= Maximum_nmax; m_prime ++ ){
            m_prime_factorial = factorial(m_prime);
            for(n =0; n <= Maximum_nmax; n++){
                n_factorial = factorial(n);
                for(m = 0; m <= Maximum_nmax; m++) {
                    m_factorial = factorial(m);

                    for(l_prime = 0; l_prime <= n_prime ; l_prime ++ ){
                        l_prime_factorial = factorial(l_prime);
                        n_prime_minus_l_prime_factorial = factorial( n_prime - l_prime );
                        for(p_prime = 0; p_prime <= m_prime; p_prime ++ ){
                            p_prime_factorial = factorial(p_prime);
                            m_prime_minus_p_prime_factorial = factorial(m_prime - p_prime);
                            for(l = 0; l<= n ; l++){
                                l_factorial = factorial(l);
                                n_minus_l_factorial = factorial(n - l);
                                for(p = 0; p<= m; p++){


                                    p_factorial = factorial(p);
                                    m_minus_p_factorial = factorial( m - p);

                                    if( Duschrinsky_rotation_Overlap[l_prime][p_prime][l][p] != 0 ){
                                        Coeff1_x = pow(A_prime, n_prime - l_prime) *
                                                sqrt( n_prime_factorial / l_prime_factorial   ) * 1/  n_prime_minus_l_prime_factorial ;
                                        Coeff1_y = pow(B_prime, m_prime - p_prime) *
                                                sqrt( m_prime_factorial / p_prime_factorial  ) * 1 /  m_prime_minus_p_prime_factorial ;

                                        Coeff2_x = pow( C_prime, n - l ) *
                                                sqrt(n_factorial /  l_factorial  )  / n_minus_l_factorial;

                                        Coeff2_y = pow(D_prime , m - p) *
                                                sqrt(m_factorial /  p_factorial  ) / m_minus_p_factorial;

                                        Coeff_prod = Coeff1_x * Coeff1_y * Coeff2_x * Coeff2_y ;
                                        if(Coeff_prod != 0){
                                            shifted_Duschrinsky_rotation_Overlap[n_prime][m_prime][n][m] =
                                                    shifted_Duschrinsky_rotation_Overlap[n_prime][m_prime][n][m] +
                                                    Coeff_prod *
                                                    Duschrinsky_rotation_Overlap[l_prime][p_prime][l][p];

                                        }
                                    }

                                }

                            }

                        }

                    }

                }
            }
        }
    }


    for(n_prime = 0; n_prime <= Maximum_nmax ; n_prime ++){
        for(m_prime = 0 ;  m_prime <= Maximum_nmax; m_prime ++ ){
                for(n =0; n <= Maximum_nmax; n++){
                    for(m = 0; m<= Maximum_nmax; m++){
                        shifted_Duschrinsky_rotation_Overlap [n_prime][m_prime][n][m] =
                                shifted_Duschrinsky_rotation_Overlap[n_prime][m_prime][n][m] * Overlap_coeff;

                    }
            }
        }
    }



}
