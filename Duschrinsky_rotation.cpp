//
// Created by phyzch on 5/8/21.
//
#include "util.h"
#include "system.h"

void detector::compute_Duschrinsky_rotation_overlap(double rotation_angle){
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


    Duschrsinky_rotation_Overlap.resize(Maximum_nmax);
    for(n_prime=0;n_prime<Maximum_nmax;n_prime++){
        Duschrsinky_rotation_Overlap[ n_prime ].resize(Maximum_nmax);
        for(m_prime=0; m_prime<Maximum_nmax; m_prime++ ){
            Duschrsinky_rotation_Overlap[ n_prime ][ m_prime ].resize(Maximum_nmax);
            for(n=0;n<Maximum_nmax;n++){
                Duschrsinky_rotation_Overlap[ n_prime ][ m_prime ][ n ].resize(Maximum_nmax);
                for(m =0; m <Maximum_nmax; m++ ){
                    Duschrsinky_rotation_Overlap[ n_prime ][ m_prime ][ n ][ m ] = 0 ;
                }
            }
        }
    }

    // same ground state
    Duschrsinky_rotation_Overlap[0][0][0][0] = 1;
    // compute [0 , m', n, m]
    for(m_prime=1; m_prime< Maximum_nmax; m_prime++){
        // compute Duschrinsky_rotation_Overlap[0][m']:
        Nonzero_element_upper_cutoff = min( m_prime - 1, Maximum_nmax - 1 );
        // n, m is index for Duschrinsky_rotation_Overlap[0][m'-1]
        // b_{y}'^{+} operator
        for( n = 0; n <= Nonzero_element_upper_cutoff ; n ++){
            for(m = 0; m <= Nonzero_element_upper_cutoff ; m ++){
                // we only have to compute contribution from nonzero coefficient.
                if(Duschrsinky_rotation_Overlap[n_prime - 1][m_prime][n][m] == 0){
                    continue ;
                }

                // cos(\theta) b_{y}^{+} operate upon |0, m'-1> .
                index_x = n ;
                index_y = m + 1;
                if(  (Maximum_nmax > index_x) and (index_x >= 0) and (Maximum_nmax > index_y) and (index_y >= 0 )   ){
                    Duschrsinky_rotation_Overlap[0][m_prime][index_x][index_y] =
                            Duschrsinky_rotation_Overlap[0][m_prime][index_x][index_y] +
                            sqrt( double(index_y) ) / sqrt( double(m_prime) ) * cos(rotation_angle) *
                            Duschrsinky_rotation_Overlap[0][m_prime-1][n][m];
                }

                // bx^{+} operator
                index_x = n + 1;
                index_y = m;
                if( (Maximum_nmax > index_x) and (index_x >= 0) and (Maximum_nmax > index_y) and (index_y >= 0 )   ){
                    Duschrsinky_rotation_Overlap[0][m_prime][index_x][index_y] =
                            Duschrsinky_rotation_Overlap[0][m_prime][index_x][index_y] +
                            sqrt( double(index_x) ) / sqrt( double(m_prime) ) * (- r_plus * sin(rotation_angle) ) *
                            Duschrsinky_rotation_Overlap[0][m_prime - 1][n][m];
                }

                // bx operator
                index_x = n - 1;
                index_y = m;
                if( (Maximum_nmax > index_x) and (index_x >= 0) and (Maximum_nmax > index_y) and (index_y >= 0 )   ){
                    Duschrsinky_rotation_Overlap[0][m_prime][index_x][index_y] =
                            Duschrsinky_rotation_Overlap[0][m_prime][index_x][index_y] +
                            sqrt( double(n) ) / sqrt( double(m_prime) ) * (- r_minus * sin(rotation_angle)) *
                            Duschrsinky_rotation_Overlap[0][m_prime - 1][n][m];
                }

            }
        }

    }

    // compute [n', m', n, m]  n' = 1, 2, 3, etc. m'= 0, 1, 2, 3, etc.
    for(n_prime = 1; n_prime < Maximum_nmax; n_prime ++ ){
        for(m_prime = 0; m_prime < Maximum_nmax; m_prime ++ ){
            Nonzero_element_upper_cutoff = min( n_prime + m_prime -1, Maximum_nmax - 1);
            // n, m is index for Duschrinsky_rotation_Overlap[n' -1][m']
            // b_{x}^{+}' operator
            for(n = 0 ; n <= Nonzero_element_upper_cutoff; n++ ){
                for(m =0; m <= Nonzero_element_upper_cutoff; m++ ){
                    // we only have to compute contribution from nonzero coefficient.
                    if(Duschrsinky_rotation_Overlap[n_prime - 1][m_prime][n][m] == 0){
                        continue ;
                    }

                    // b_{x}^{+} operator
                    index_x = n + 1;
                    index_y = m;
                    if(  (Maximum_nmax > index_x) and (index_x >= 0) and (Maximum_nmax > index_y) and (index_y >= 0 )   ){
                        Duschrsinky_rotation_Overlap[n_prime][m_prime][index_x][index_y] =
                                Duschrsinky_rotation_Overlap[n_prime][m_prime][index_x][index_y] +
                                sqrt( double(index_x) ) / sqrt(  double(n_prime)  ) * cos(rotation_angle) *
                                Duschrsinky_rotation_Overlap[n_prime-1][m_prime][n][m];
                    }

                    // b_{y}^{+} operator
                    index_x = n;
                    index_y = m + 1;
                    if(  (Maximum_nmax > index_x) and (index_x >= 0) and (Maximum_nmax > index_y) and (index_y >= 0 )   ){
                        Duschrsinky_rotation_Overlap[n_prime][m_prime][index_x][index_y] =
                                Duschrsinky_rotation_Overlap[n_prime][m_prime][index_x][index_y] +
                                sqrt(  double(index_y)  ) / sqrt( double(n_prime) ) * r_plus * sin(rotation_angle) *
                                Duschrsinky_rotation_Overlap[n_prime - 1][m_prime][n][m];
                    }

                    // b_{y} operator
                    index_x = n;
                    index_y = m - 1;
                    if(  (Maximum_nmax > index_x) and (index_x >= 0) and (Maximum_nmax > index_y) and (index_y >= 0 )   ){
                        Duschrsinky_rotation_Overlap[n_prime][m_prime][index_x][index_y] =
                                Duschrsinky_rotation_Overlap[n_prime][m_prime][index_x] [index_y] +
                                sqrt( double(m) ) / sqrt(double(n_prime)) * (- r_minus) * sin(rotation_angle) *
                                Duschrsinky_rotation_Overlap[n_prime - 1] [m_prime][n][m];
                    }

                }
            }

        }

    }


}

void detector:: shifted_Duschrinsky_rotation_overlap(){
    
}
