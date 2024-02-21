#include <stdio.h>
#include <math.h>

/* prototype the functions */
void mat3x3_inv(double minv[][3], const double m[][3]);
void mat3x3_vec_mult(double b[], const double A[][3], const double x[]);
void ikine(double theta[], const double x[], const double y[], double phi[], const double L[]);
void jacob(double m[][3], const double theta[], const double L[]);

int main(int argc, char *argv[])
{
    /* introduce yourself */
    printf("Hello World! How are you today?\n\n");

    /* type the variables */
    double initpose[3];
    double finalpose[3];
    double linklength[3];
    double time = 2;
    const unsigned int N = 10;
    double J[3][3];
    double Ji[3][3];
    double endpoint_velocity[3];
    double omega[3];
    double theta[3];
    double x[N+2];
    double y[N+2];
    double phi[N+2];
    double dx;
    double dy;
    double dphi;

    /* define initial pose */
    initpose[0] = 3;
    initpose[1] = 3;
    initpose[2] = 0;

    /* define final pose */
    finalpose[0] = 5;
    finalpose[1] = 5;
    finalpose[2] = 0;

    /* define link lengths of serial manipulator */
    linklength[0] = 5;
    linklength[1] = 4;
    linklength[2] = 3;

    /* define average velocity goals of the endpoint */
    endpoint_velocity[0] = (finalpose[0] - initpose[0])/time;
    endpoint_velocity[1] = (finalpose[1] - initpose[1])/time;
    endpoint_velocity[2] = (finalpose[2] - initpose[2])/time;

    /* sample the endpoints on the trajectory */
    dx = (finalpose[0] - initpose[0])/(N+1);
    dy = (finalpose[1] - initpose[1])/(N+1);
    dphi = (finalpose[2] - initpose[2])/(N+1);

    x[0] = initpose[0];
    y[0] = initpose[1];
    phi[0] = initpose[2];

    x[N+1] = finalpose[0];
    y[N+1] = finalpose[1];
    phi[N+1] = finalpose[2];

    for (unsigned int i=1; i<N+1; i++)
    {
        x[i] = x[i-1] + dx;
        y[i] = y[i-1] + dy;
        phi[i] = phi[i-1] + dphi;
    }

    printf("the trajectory of sampled endpoint locations are\n");
    printf("================================================\n");
    for (unsigned int j=0; j<N+2; j++)
    {
        printf("(x, y, phi) = (%f, %f, %f)\n", x[j],y[j],phi[j]);
    }

    /* compute the inverse kinematics and joint speeds */
    printf("\n");
    printf("the corresponding joint-angles and joint-speeds are\n");
    printf("===================================================\n");
    for (unsigned int k=0; k<N+2; k++)
    {
        ikine(theta, x+k, y+k, phi+k, linklength);
        jacob(J, theta, linklength);
        mat3x3_inv(Ji, J);
        mat3x3_vec_mult(omega,Ji,endpoint_velocity);
        printf("(theta1, theta2, theta3) \t (omega1, omega2, omega3)");
        printf(" = (%f, %f, %f) \t (%f, %f, %f)\n", theta[0], theta[1], theta[2], omega[0], omega[1], omega[2]);
    }

    return 0;
}

/* inverse kinematics of 3-link planar serial manipulator */
void ikine(double theta[], const double x[], const double y[], double phi[], const double L[])
{
    double xb = x[0] - (L[2] * cos(phi[0]));
    double yb = y[0] - (L[2] * sin(phi[0]));

    double num = pow((L[0]+L[1]),2) - (pow(xb,2)+pow(yb,2));
    double den = (pow(xb,2)+pow(yb,2)) - pow((L[0]-L[1]),2);
    double sqrtval = sqrt(num/den);

    theta[1] = 2*atan(sqrtval);
    theta[0] = atan2(yb,xb) - atan2(L[1]*sin(theta[1]), L[0]+L[1]*cos(theta[1]));
    theta[2] = phi[0] - theta[0] - theta[1];
}

/* jacobian matrix */
void jacob(double m[][3], const double theta[], const double L[])
{
    double s1 = sin(theta[0]);
    double c1 = cos(theta[0]);
    double s2 = sin(theta[0]+theta[1]);
    double c2 = cos(theta[0]+theta[1]);
    double s3 = sin(theta[0]+theta[1]+theta[2]);
    double c3 = cos(theta[0]+theta[1]+theta[2]);

    /* row one */
    m[0][0] = -L[0]*s1 - L[1]*s2 - L[2]*s3;
    m[0][1] = -L[1]*s2 - L[2]*s3;
    m[0][2] = -L[2]*s3;

    /* row two */
    m[1][0] = L[0]*c1 + L[1]*c2 + L[2]*c3;
    m[1][1] = L[1]*c2 + L[2]*c3;
    m[1][2] = L[2]*c3;

    /* row three */
    m[2][0] = 1;
    m[2][1] = 1;
    m[2][2] = 1;
}

/* Direct inverse of a 3x3 matrix m */
void mat3x3_inv(double minv[][3], const double m[][3])
{
    double det;
    double invdet;

    det =   m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    invdet = 1.0 / det;

    minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
    minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;

    minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
    minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;

    minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
    minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
    minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}

/* Explicit 3x3 matrix * vector b = A x where A is 3x3 */
void mat3x3_vec_mult(double b[], const double A[][3], const double x[])
{
    const int N = 3;
    int idx, jdx;

    for (idx = 0; idx < N; idx++)
    {
        b[idx] = 0.0;

        for (jdx = 0; jdx < N; jdx++)
        {
            b[idx] += A[idx][jdx] * x[jdx];
        }
    }
}
