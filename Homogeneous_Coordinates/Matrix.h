#pragma once
#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <array>
#include <type_traits>



template <typename T, size_t N, size_t M>
class Matrix {
private:
    std::array<T, N* M> m_matrix;
public:

    Matrix() = default;
    ~Matrix() = default;

    T operator[] (size_t i) const;
    T& operator[] (size_t i);

    //operation with scalars
    template <typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
    auto operator+ (const T2 add);

    template <typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
    auto operator- (const T2 add);

    template <typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
    auto operator* (const T2 add);

    template <typename T2, typename = std::enable_if_t<std::is_arithmetic_v<T2>>>
    auto operator/ (const T2 add);

    //template <typename T2>
    //friend auto operator= (const T2 add);

    //unitary operation 
    Matrix<T, N, M>& operator++ (); //prefix

    Matrix<T, N, M> operator++ (int); //postfix

    Matrix<T, N, M>& operator-- (); //prefix

    Matrix<T, N, M> operator-- (int); //postfix

    //binary operation between matrix
    template <typename T2>
    auto operator+ (const Matrix<T2, N, M>& add);

    template <typename T2>
    auto operator- (const Matrix<T2, N, M>& sub);

    template<typename T2, size_t M2>
    auto operator* (const Matrix<T2, M, M2>& mult);

    //boolean
    template<typename T2, size_t N2, size_t M2>
    int operator ==(const Matrix<T2, N2, M2>& i_matrix);


    Matrix<T, N, M>& operator= (const Matrix<T, N, M>& copy);


    //mathods

    size_t Height() const;
    size_t Width() const;

    template <typename = std::enable_if_t<N == M>>
    T  Determinant() const;

    Matrix<T, M, N> Transpose();
};



template <typename T, size_t N, size_t M>
size_t Matrix<T, N, M>::Height() const
{
    return N;
}
template <typename T, size_t N, size_t M>
size_t Matrix<T, N, M>::Width() const
{
    return M;
}

template <typename T, size_t N, size_t M>
template <typename>
T  Matrix<T, N, M>::Determinant() const
{
    if constexpr (N == 1)
    {
        std::cout << "Determinant for N=M = 1\n";
        //return (m_matrix[0] < 0 ? (-m_matrix[0]) : m_matrix[0]);
        return  m_matrix[0];
    }
    else if (N == 2)
    {
        std::cout << "Determinant for N=M = 2\n";
        return m_matrix[0] * m_matrix[3] - m_matrix[1] * m_matrix[2];
    }
    else if (N == 3)
    {
        std::cout << "Determinant for N=M = 3\n";
        return (m_matrix[0] * m_matrix[4] * m_matrix[8] + 
            m_matrix[1] * m_matrix[5] * m_matrix[6] +
            m_matrix[2] * m_matrix[3] * m_matrix[7]) -
            (m_matrix[0] * m_matrix[5] * m_matrix[7] +
            m_matrix[1] * m_matrix[3] * m_matrix[8] +
            m_matrix[2] * m_matrix[4] * m_matrix[6]);
    }
    else
    {
        std::cout << "Determinant for N=M >= 4 (Has no realisation)\n";
        return  static_cast<T>(0);
    }
}

//operation with scalars
template <typename T, size_t N, size_t M>
template <typename T2, typename>
auto Matrix<T, N, M>::operator+ (const T2 number)
{
    using resultType = decltype(std::declval<T>() + std::declval<T2>());
    Matrix<resultType, N, M> result;
    for (int i = 0; i < N * M; ++i)
    {
        result[i] = this->m_matrix[i] + number;
    }

    return result;
}

template <typename T, size_t N, size_t M>
template <typename T2, typename>
auto Matrix<T, N, M>::operator- (const T2 number)
{
    using resultType = decltype(std::declval<T>() - std::declval<T2>());
    Matrix<resultType, N, M> result;
    for (int i = 0; i < N * M; ++i)
    {
        result[i] = this->m_matrix[i] - number;
    }

    return result;
}

template <typename T, size_t N, size_t M>
template <typename T2, typename>
auto Matrix<T, N, M>::operator* (const T2 number)
{
    using resultType = decltype(std::declval<T>()* std::declval<T2>());;
    Matrix<resultType, N, M> result;
    for (int i = 0; i < N * M; ++i)
    {
        result[i] = this->m_matrix[i] * number;
    }

    return result;
}

template <typename T, size_t N, size_t M>
template <typename T2, typename>
auto Matrix<T, N, M>::operator/ (const T2 number)
{
    using resultType = decltype(std::declval<T>() / std::declval<T2>());
    Matrix<resultType, N, M> result;
    for (int i = 0; i < N * M; ++i)
    {
        result[i] = this->m_matrix[i] / number;
    }

    return result;
}

//template <typename T2>
//auto operator= (const T2 number);


//unitary operation 
template <typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator++ () ////prefix
{
    for (int i = 0; i < N * M; ++i)
    {
        this->m_matrix[i]++;
    }

    return *this;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M> Matrix<T, N, M>::operator++ (int) //postfix
{
    Matrix<T, N, M> temp = *this;
    for (int i = 0; i < N * M; ++i)
    {
        this->m_matrix[i]++;
    }

    return temp;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator-- () //prefix
{
    for (int i = 0; i < N * M; ++i)
    {
        this->m_matrix[i]--;
    }

    return *this;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M> Matrix<T, N, M>::operator-- (int) //postfix
{
    Matrix<T, N, M> temp = *this;
    for (int i = 0; i < N * M; ++i)
    {
        this->m_matrix[i]--;
    }

    return temp;
}

//binary operation between matrix
template <typename T, size_t N, size_t M>
template <typename T2>
auto Matrix<T, N, M>::operator+ (const Matrix<T2, N, M>& add)
{
    using resultType = decltype(std::declval<T>() + std::declval<T2>());
    Matrix<resultType, N, M> result;

    for (int i = 0; i < N * M; ++i)
        result[i] = this->m_matrix[i] + add[i];

    return result;
}

template <typename T, size_t N, size_t M>
template <typename T2>
auto Matrix<T, N, M>::operator- (const Matrix<T2, N, M>& sub)
{
    using resultType = decltype(std::declval<T>() - std::declval<T2>());
    Matrix<resultType, N, M> result;

    for (int i = 0; i < N * M; ++i)
        result[i] = this->m_matrix[i] - sub[i];

    return result;
}


template<typename T, size_t N, size_t M>
template<typename T2, size_t M2>
auto Matrix<T, N, M>::operator*(const Matrix<T2, M, M2>& mult)
{
    using resultType = decltype(std::declval<T>()* std::declval<T2>());;

    Matrix<resultType, N, M2> result;

    resultType temp;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < M2; j++)
        {
            temp = 0;
            for (size_t k = 0; k < M; k++)
            {
                T temp1 = this->m_matrix[i * M + k];
                T2 temp2 = mult[k * M2 + j];
                temp += this->m_matrix[i * M + k] * mult[k * M2 + j];
            }
            result[i * M2 + j] = temp;
        }
    }
    return result;
}

template <typename T, size_t N, size_t M>
Matrix<T, N, M>& Matrix<T, N, M>::operator= (const Matrix<T, N, M>& copy)
{
    this->m_matrix = copy.m_matrix;
    return *this;
}

template <typename T, size_t N, size_t M>
T  Matrix<T, N, M>::operator[](size_t i) const
{
    return m_matrix[i];
}

template <typename T, size_t N, size_t M>
T& Matrix<T, N, M>::operator[](size_t i)
{
    return m_matrix[i];
}



template <typename T, size_t N, size_t M>
std::ostream& operator<<(std::ostream& os, const Matrix<T, N, M>& out)
{
    if constexpr (N == 1)
    {
        for (size_t j = 0; j < M; j++)
            os << out[j] << ' ';
        os << '\n';
    }
    else if (M == 1)
    {
        for (size_t i = 0; i < N; i++)
            os << out[i] << '\n';
    }
    else
    {
      //  os << "Output for else\n";
        for (size_t i = 0; i < N; i++)
        {
            for (size_t j = 0; j < M; j++)
                os << out[i * M + j] << ' ';
            os << '\n';
        }
    }
    

    return os;
}

template <typename T, size_t N, size_t M>
template<typename T2, size_t N2, size_t M2>
int Matrix<T, N, M>::operator ==(const Matrix<T2, N2, M2>& i_matrix)
{
    if constexpr (N != N2 || M != M2)
    {
        return 0;
    }
    else
    {
        for (size_t i = 0; i < N*M; i++)
        {
            if (this->m_matrix[i] != i_matrix[i])
                return 0;

        }
        return 1;
    }
}

template <typename T, size_t N, size_t M>
Matrix<T, M, N> Matrix<T, N, M>::Transpose()
{
    if constexpr (N == 1 || M == 1)
    {
        Matrix<T, M, N> result;
        for (size_t i = 0; i < M*N; i++)
            result[i] = this->m_matrix[i];
        
        return result;
    }
    else
    {
        Matrix<T, M, N> result;
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < M; j++)
            {
                size_t oldInd = i * M + j;
                size_t newInd = j * N + i;
                result[j * N + i] = this->m_matrix[i * M + j];
            }
                

        return result;
    }
}

