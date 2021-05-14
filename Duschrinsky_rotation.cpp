//
// Created by phyzch on 5/8/21.
//
#include "util.h"
#include "system.h"

void detector::compute_ground_state_overlap(){
    int i, j;
    double omega_y = mfreq[0][2];
    double omega_x = mfreq[0][1];

    double r_plus = double(1)/2 * (sqrt(omega_y / omega_x)  + sqrt(omega_x / omega_y) );
    double r_minus = double(1)/2 * (sqrt(omega_y / omega_x) - sqrt(omega_x / omega_y) );

    double A = cos(rotation_angle);
    double B = r_plus * sin(rotation_angle);
    double C = - r_minus * sin(rotation_angle);

    int Maximum_nmax = max( nmax[0][1], nmax[0][2] );
    double norm = 0;
    double overlap_cutoff = 1e-2;

    for(i=0; i< Maximum_nmax;i++){
        vector<double> v1;
        for(j=0; j< Maximum_nmax; j++){
            v1.push_back(0);
        }
        ground_state_overlap.push_back(v1);
    }

    ground_state_overlap[0][0] = 1;
#   // finish [0][^] , [^][0] and [1,^] [^,1] terms
    for(i=2;i<Maximum_nmax - 2 ; i = i + 2){
        if(ground_state_overlap[i-2][0]!=0){
            ground_state_overlap[i][0] = sqrt(double(i-1)/ i ) * B * C / ( pow(A,2) + pow(B,2) ) * ground_state_overlap[i-2][0];
            ground_state_overlap[i-1][1] = - sqrt(i-1) * A * C * (pow(A,2) + pow(B,2) ) * ground_state_overlap[i-2][0];

            if( abs(ground_state_overlap[i][0]) < overlap_cutoff){
                ground_state_overlap[i][0] = 0;
            }
            if( abs(ground_state_overlap[i-1][1]) < overlap_cutoff){
                ground_state_overlap[i-1][1] = 0;
            }

        }

    }
    for(i=2; i<Maximum_nmax - 2; i= i+2){
        if(ground_state_overlap[0][i-2]!=0){
            ground_state_overlap[0][i] = - sqrt( double(i-1)/i ) * B * C / (pow(A,2) + pow(B,2) ) * ground_state_overlap[0][i-2];
            ground_state_overlap[1][i-1] = -sqrt(i-1) * A * C / (pow(A,2) + pow(B,2)) * ground_state_overlap[0][i-2];

            if( abs(ground_state_overlap[0][i]) < overlap_cutoff){
                ground_state_overlap[0][i] = 0;
            }
            if( abs(ground_state_overlap[1][i-1]) < overlap_cutoff){
                ground_state_overlap[1][i-1] = 0;
            }
        }

    }

    for(i=1; i< Maximum_nmax - 1; i++){
        for(j=i;j+2 < Maximum_nmax; j=j+2){
            ground_state_overlap[i+1][j+1] = - (C * sqrt(j+1) * ground_state_overlap[i][j] + B * sqrt(j+2) * ground_state_overlap[i][j+2]) / (A * sqrt(i+1));
            if( abs(ground_state_overlap[i+1][j+1]) < overlap_cutoff){
                ground_state_overlap[i+1][j+1] = 0;
            }

        }
        for(j=i;j+2 < Maximum_nmax; j = j+2){
            ground_state_overlap[j+1][i+1] = (B * sqrt(j+2) * ground_state_overlap[j+2][i] - C * sqrt(j+2) * ground_state_overlap[j][i]) / (A * sqrt(i+1)) ;
            if( abs(ground_state_overlap[j+1][i+1]) < overlap_cutoff){
                ground_state_overlap[j+1][i+1] = 0;
            }

        }
    }

    // normalize overlap.
    for(i=0;i<Maximum_nmax;i++){
        for(j=0;j<Maximum_nmax; j++){
            norm = norm + pow(ground_state_overlap[i][j] , 2);
        }
    }
    norm = sqrt(norm);
    for(i=0;i<Maximum_nmax;i++){
        for(j=0;j<Maximum_nmax; j++){
            ground_state_overlap[i][j] = ground_state_overlap[i][j] / norm ;
        }
    }

    for(i=0;i<Maximum_nmax;i++){
        for(j=0;j<Maximum_nmax;j++){
            if(ground_state_overlap[i][j] != 0){
                nonzero_ground_state_overlap_coeff.push_back(ground_state_overlap[i][j]);
                vector<int> v {i , j};
                nonzero_ground_state_overlap_index.push_back(v);
            }
        }
    }

}

void detector::compute_Duschrinsky_rotation_overlap(){
    // mode 1 and mode 2 is vibrational coordinate we perform Duschrinsky rotation
    int i, j;
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

    int max_index = 0;
    int Len = nonzero_ground_state_overlap_index.size();
    for(i=0;i<Len; i ++ ){
        for(j=0;j<2;j++){
            if (max_index < nonzero_ground_state_overlap_index[i][j]){
                max_index = nonzero_ground_state_overlap_index[i][j];
            }
        }
    }

    Duschrinsky_rotation_Overlap.resize(Maximum_nmax + 1);
    for(n_prime=0;n_prime <= Maximum_nmax;n_prime++){
        Duschrinsky_rotation_Overlap[ n_prime ].resize(Maximum_nmax + 1);
        for(m_prime=0; m_prime <= Maximum_nmax; m_prime++ ){
            Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ].resize(2 * Maximum_nmax + 1 + max_index  );
            for(n=0;n <= 2 * Maximum_nmax + max_index ;n++){
                Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ][ n ].resize(2 * Maximum_nmax + 1 + max_index );
                for(m =0; m <= 2 * Maximum_nmax + max_index ; m++ ){
                    Duschrinsky_rotation_Overlap[ n_prime ][ m_prime ][ n ][ m ] = 0 ;
                }
            }
        }
    }


    // ground state for two oscillator is not the same. Need to fix bug here. Although [0,0] is major contribution
    for(i=0;i<Len;i++){
        index_x = nonzero_ground_state_overlap_index[i][0];
        index_y = nonzero_ground_state_overlap_index[i][1];
        Duschrinsky_rotation_Overlap[0][0][index_x][index_y] = nonzero_ground_state_overlap_coeff[i];
    }
    // compute [0 , m', n, m]
    for(m_prime=1; m_prime<= Maximum_nmax; m_prime++){
        // compute Duschrinsky_rotation_Overlap[0][m']:
        Nonzero_element_upper_cutoff =  m_prime - 1 + max_index ;
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
            Nonzero_element_upper_cutoff = n_prime + m_prime -1 + max_index;
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
