#include "Matrix.h"

#include <iostream>

using fPoint = float;

template class Matrix<fPoint, 3, 3>;


inline Matrix<fPoint, 3, 1> transform(const Matrix<fPoint, 3, 1> &i_coordToTransf,  Matrix<fPoint, 3,3> &i_trMatrix)
{
    return i_trMatrix * i_coordToTransf;
}

//normalise
inline Matrix<fPoint, 3, 1> normalise(const Matrix<fPoint, 3, 1>&& i_coordToTransf)
{
    Matrix<fPoint, 3, 1> normalise = i_coordToTransf;
    for (size_t i = 0; i < 3; i++)
    {
        normalise[i] = i_coordToTransf[i]/i_coordToTransf[2];
    }
    return normalise;
}


/*
* 
^                           ^                 p3
|                           |             ~~~+
|                           |    p1   ~~~~    \
|                           |     +~~~         \
+--------+           ==>    |    /              \
|        |                  |   /                \
|        |                  |  +~~~~~~~~~         \
|        |                  | p0         ~~~~~~~~~~+ p2
+--------+------->          *--------------------------->

*/
// transformation matrix to transform singular 2D cube((0;0),(0;1),(1;0),(1;1)) into quadrilateral(p0,p1,p2,p3) 
Matrix<fPoint, 3, 3> foundTransformationMatrixLecture(Matrix<fPoint, 3, 1>& p0,
                                                Matrix<fPoint, 3, 1>& p1,
                                                Matrix<fPoint, 3, 1>& p2,
                                                Matrix<fPoint, 3, 1>& p3)
{
    fPoint d, a, b;
    fPoint x0 = p0[0], y0 = p0[1], x1 = p1[0], y1 = p1[1],
        x2 = p2[0], y2 = p2[1], x3 = p3[0], y3 = p3[1];

    d = 1 / ((x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3));
    a = ((y1 - y3) * (x2 - 2 * x0 + 1) + (y0 - y2) * (x1 - x3)) * d;
    b = ((x2 - x3) * (-y0 + y1 + y2 - y3) - (y2 - y3) * (-2 * x0 + x1 + x2 - x3 + 1)) * d;

    Matrix<fPoint, 3, 3> result;

    result[0] = a * x2 - x0 + x2;   //A
    result[1] = b * x1 - x0 + x1;   //B
    result[2] = x0;                 //C
    result[3] = a * y2 - y0 + y2;   //D
    result[4] = b * y1 - y0 + y1;   //E
    result[5] = y0;                 //F
    result[6] = a;                  //a
    result[7] = b;                  //b
    result[8] = 1;                  //c

    return result;
}



// transformation matrix to transform singular 2D cube((0;0),(0;1),(1;0),(1;1)) into quadrilateral(p0,p1,p2,p3) 
Matrix<fPoint, 3, 3> foundTransformationMatrixMy(const Matrix<fPoint, 3, 1>& p0,
                                                const Matrix<fPoint, 3, 1>& p1,
                                                const Matrix<fPoint, 3, 1>& p2,
                                                const Matrix<fPoint, 3, 1>& p3)
{
    fPoint d, a, b;
    fPoint x0 = p0[0], y0 = p0[1], x1 = p1[0], y1 = p1[1], x2 = p2[0], y2 = p2[1], x3 = p3[0], y3 = p3[1];

    d = 1 / (x1 * y2 - x1 * y3 - x2 * y1 + x2 * y3 + x3 * y1 - x3 * y2);
    a = (-x0 * y1 + x0 * y3 + x1 * y0 - x1 * y2 + x2 * y1 - x2 * y3 - x3 * y0 + x3 * y2) * d;
    b = (x0 * y2 - x0 * y3 - x1 * y2 + x1 * y3 - x2 * y0 + x2 * y1 + x3 * y0 - x3 * y1) * d;;

    Matrix<fPoint, 3, 3> result;
    result[0] = a * x2 - x0 + x2;   //A
    result[1] = b * x1 - x0 + x1;   //B
    result[2] = x0;                 //C
    result[3] = a * y2 - y0 + y2;   //D
    result[4] = b * y1 - y0 + y1;   //E
    result[5] = y0;                 //F
    result[6] = a;                  //a
    result[7] = b;                  //b
    result[8] = 1;                  //c

    return result;
}


/*
                                 
^      + p1                       ^
|     / \                         |
|    /   \                        +
|   /     \             ==>       | \
|  +~~~~   \                      |    \
| p0    ~~~~+ p2                  |       \
*----------------->               +---------+-->

*/
// transformation matrix to transform any triangle(p0,p1,p2) into ((0;0),(1;0),(0;1)) 
Matrix<fPoint, 3, 3> foundTransformationMatrix(const Matrix<fPoint, 3, 1>& p0,
                                                const Matrix<fPoint, 3, 1>& p1,
                                                const  Matrix<fPoint, 3, 1>& p2)
{
    fPoint d, x0 = p0[0], y0 = p0[1], x1 = p1[0], y1 = p1[1], x2 = p2[0], y2 = p2[1];

    Matrix<fPoint, 3, 3> result;
    d = 1 / (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1));
    //d = 1 / (x0 * y1 - x0 * y2 - x1 * y0 + x1 * y2 + x2 * y0 - x2 * y1);
    result[0] = (y2 - y0) * d;              //A
    result[1] = (x0 - x2) * d;              //B
    result[2] = (x2 * y0 - x0 * y2) * d;    //C
    result[3] = (y0 - y1) * d;              //D
    result[4] = (x1 - x0) * d;              //E
    result[5] = (x0 * y1 - x1 * y0) * d;    //F
    result[6] = 0;                          //a
    result[7] = 0;                          //b
    result[8] = 1;                          //c

    return result;
}



int main()
{

    //test to transform singular 2D cube((0;0),(0;1),(1;0),(1;1)) into quadrilateral(p0,p1,p2,p3)  
    Matrix<fPoint, 3, 1> p0; p0[0] = 1; p0[1] = 1; p0[2] = 1;
    Matrix<fPoint, 3, 1> p1; p1[0] = 0; p1[1] = 3; p1[2] = 1;
    Matrix<fPoint, 3, 1> p2; p2[0] = 3; p2[1] = 1; p2[2] = 1;
    Matrix<fPoint, 3, 1> p3; p3[0] = 4; p3[1] = 3; p3[2] = 1;

    Matrix<fPoint, 3, 1> p_0; p_0[0] = 0; p_0[1] = 0; p_0[2] = 1;
    Matrix<fPoint, 3, 1> p_1; p_1[0] = 0; p_1[1] = 1; p_1[2] = 1;
    Matrix<fPoint, 3, 1> p_2; p_2[0] = 1; p_2[1] = 0; p_2[2] = 1;
    Matrix<fPoint, 3, 1> p_3; p_3[0] = 1; p_3[1] = 1; p_3[2] = 1;

    auto trMatrLec = foundTransformationMatrixLecture(p0, p1, p2, p3);
    auto trMatrMy = foundTransformationMatrixMy(p0, p1, p2, p3);

    if (trMatrLec == trMatrMy)
    {
        std::cout << "Matrixes are equal. Transformation matrix:\n" << trMatrLec << '\n';

        std::cout << "trMatrLec*p0 = " << normalise((trMatrLec * p_0)).Transpose()
                << "trMatrLec*p1 = " << normalise((trMatrLec * p_1)).Transpose()
                << "trMatrLec*p2 = " << normalise((trMatrLec * p_2)).Transpose()
                << "trMatrLec*p3 = " << normalise((trMatrLec * p_3)).Transpose()
                << '\n';
    }
    else
    {
        std::cout << "foundTransformationMatrixLecture:\n" << trMatrLec
            << "trMatrLec*p0 = " << trMatrLec * p_0
            << "trMatrLec*p1 = " << trMatrLec * p_1
            << "trMatrLec*p2 = " << trMatrLec * p_2
            << "trMatrLec*p3 = " << trMatrLec * p_3
            << '\n';

        std::cout << "foundTransformationMatrixMy:\n" << trMatrLec
            << "trMatrMy*p0 = " << trMatrMy * p_0
            << "trMatrMy*p1 = " << trMatrMy * p_1
            << "trMatrMy*p2 = " << trMatrMy * p_2
            << "trMatrMy*p3 = " << trMatrMy * p_3;

    }


    
    //test to transform triangle((1;2),(5;1),(4;5)) into ((0;0),(1;0),(0;1))
    p0[0] = 1; p0[1] = 2; p0[2] = 1;
    p1[0] = 5; p1[1] = 1; p1[2] = 1;
    p2[0] = 4; p2[1] = 5; p2[2] = 1;

    p_0[0] = 0; p_0[1] = 0; p_0[2] = 1;
    p_1[0] = 1; p_1[1] = 0; p_1[2] = 1;
    p_2[0] = 0; p_2[1] = 1; p_2[2] = 1;
    trMatrMy = foundTransformationMatrix(p0, p1, p2);


    std::cout << "Transformation matrix for triangle:\n" << trMatrLec << '\n';
    std::cout << "trMatrLec*p0 = " << (trMatrMy * p0).Transpose()
            << "trMatrLec*p1 = " << (trMatrMy * p1).Transpose()
            << "trMatrLec*p2 = " << (trMatrMy * p2).Transpose()
            << '\n';





    return 0;
}