#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

using VECTOR = std::vector<double>;

using MASSIVE = std::vector<VECTOR>;

class MATRIX
{
public:
    MASSIVE data;
    int rows_;
    int columns_;
    // Конструкторы

    MATRIX() : rows_(0), columns_(0) {}

    MATRIX(int rows_, int columns_) : rows_(rows_), columns_(columns_), data(rows_, VECTOR(columns_)) {}

    MATRIX(const std::vector<VECTOR> &data) : data(data), rows_(data.size()), columns_(data[0].size()) {}

    MATRIX(const MATRIX &matr) : data(matr.data), rows_(matr.rows_), columns_(matr.columns_) {}

    ~MATRIX() {}

    // Геттеры, Сеттеры

    int getRows() const { return rows_; }
    int getColumns() const { return columns_; }
    int getDim() const { return getRows() * getColumns(); }
    double getElement(int row, int col) const { return data[row][col]; }
    VECTOR getRow(int row) const { return data[row]; }
    VECTOR getCol(int col) const
    {
        VECTOR colMatr(getRows());
        for (int i = 0; i < getRows(); ++i)
            colMatr[i] = data[i][col];
        return colMatr;
    }
    std::vector<VECTOR> getData() const { return data; }
    std::string toString() const
    {
        std::stringstream ss;
        for (const auto &row : data)
        {
            for (const auto &elem : row)
                ss << elem << " ";
            ss << "\n";
        }
        return ss.str();
    }

    void setElement(int row, int col, double val)
    {
        if (row >= 0 && row < rows_ && col >= 0 && col < columns_)
            data[row][col] = val;
        else
            std::cout << "Invaled row or column index." << std::endl;
    }
    void setRow(int row, const VECTOR &val)
    {
        if (row >= 0 && row < rows_ && val.size() == columns_)
        {
            for (int i = 0; i < columns_; ++i)
                data[row][i] = val[i];
        }
        else
            std::cout << "Invalid row index or value size." << std::endl;
    }
    void setCol(int col, const VECTOR &val)
    {
        if (col >= 0 && col < columns_ && val.size() == rows_)
        {
            for (int i = 0; i < rows_; ++i)
                data[i][col] = val[i];
        }
        else
            std::cout << "Invalid row index or value size." << std::endl;
    }

    // Перегрузка операторов

    double &operator()(int row, int col) { return data[row][col]; }      // изменение значения эл-та через оператор ()
    double operator()(int row, int col) const { return data[row][col]; } // получение значения эл-та через оператор ()
    MATRIX operator+(const MATRIX &matr) const
    { // сложение матриц
        if (rows_ != matr.rows_ || columns_ != matr.columns_)
            std::cout << "These matrices cannot be folded." << std::endl;
        MATRIX res(rows_, columns_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                res(i, j) = data[i][j] + matr(i, j);
        return res;
    }
    MATRIX operator-(const MATRIX &matr) const
    { // вычитание матриц
        if (rows_ != matr.rows_ || columns_ != matr.columns_)
            std::cout << "These matrices cannot be subtracted." << std::endl;
        MATRIX res(rows_, columns_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                res(i, j) = data[i][j] - matr(i, j);
        return res;
    }
    MATRIX operator*(const MATRIX &matr) const
    { // умножение матриц
        if (columns_ != matr.rows_)
            std::cout << "MATRIX dimensions are not compatible for multiplicztion." << std::endl;
        MATRIX res(rows_, matr.columns_);
        for (int i = 0; i < rows_; ++i)
        {
            for (int j = 0; j < matr.columns_; ++j)
            {
                double sum = 0;
                for (int k = 0; k < columns_; ++k)
                    sum += data[i][k] * matr(k, j);
                res(i, j) = sum;
            }
        }
        return res;
    }
    bool operator==(const MATRIX &matr) const
    { // равенство матриц
        if (rows_ != matr.rows_ || columns_ != matr.columns_)
            return false;
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                if (data[i][j] != matr(i, j))
                    return false;
        return true;
    }
    bool operator!=(const MATRIX &matr) const { return !(*this == matr); } // неравенство матриц
    MATRIX operator*(double digit) const
    { // умножение матрицы на число
        MATRIX res(rows_, columns_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                res(i, j) = data[i][j] * digit;
        return res;
    }
    MATRIX operator/(double digit) const
    { // деление матрицы на число
        if (digit == 0)
            std::cout << "The division cannot be performed." << std::endl;
        MATRIX res(rows_, columns_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                res(i, j) = data[i][j] / digit;
        return res;
    }

    // Работа с Матрицами

    void Resize(int newRows, int newColumns, bool preserve = false)
    { // изменение размерности
        if (preserve)
        {
            if (newRows < rows_ || newColumns < columns_)
                std::cout << "Cannot reduce size while preserving data." << std::endl;
            std::vector<VECTOR> newMatr(newRows, VECTOR(newColumns));
            for (int i = 0; i < rows_; ++i)
                for (int j = 0; j < columns_; ++j)
                    newMatr[i][j] = data[i][j];
            data = newMatr;
            rows_ = newRows;
            columns_ = newColumns;
        }
        else
        {
            data.resize(newRows, VECTOR(newColumns));
            rows_ = newRows;
            columns_ = newColumns;
        }
    }

    bool isE() const
    { // проверка: единичная ли матрица
        if (rows_ != columns_)
            return false;
        for (int i = 0; i < rows_; ++i)
        {
            for (int j = 0; j < columns_; ++j)
            {
                if (i == j && data[i][j] != 1)
                    return false;
                else if (i != j && data[i][j] != 0)
                    return false;
            }
        }
        return true;
    }
    bool isSqr() const { return rows_ == columns_; } // проверка: квадратная ли матрица

    double det() const
    { // определитель матрицы
        if (!isSqr())
            return 0;
        if (rows_ == 1)
            return data[0][0];
        if (rows_ == 2)
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        double det = 0.0;
        for (int i = 0; i < columns_; ++i)
            det += (i % 2 == 0 ? 1 : -1) * data[0][i] * MMatrix(0, i).det();
        return det;
    }
    MATRIX transp() const
    { // транспонированная матрица
        MATRIX res(columns_, rows_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                res(j, i) = data[i][j];
        return res;
    }
    MATRIX MMatrix(int row, int col) const
    { // матрица минора для эл-та по индексу
        MATRIX res(rows_ - 1, columns_ - 1);
        int ri = 0, rj = 0;
        for (int i = 0; i < rows_; ++i)
        {
            if (i != row)
            {
                for (int j = 0; j < columns_; ++j)
                {
                    if (j != col)
                        res(ri, rj++) = data[i][j];
                }
                rj = 0;
                ri++;
            }
        }
        return res;
    }
    double Minordet(int row, int col) const
    { // определитель матрицы минора для эл-та по индексу
        return MMatrix(row, col).det();
    }
    MATRIX MatrixMinors() const
    { // матрица минора для соответствующих эл-там и их индексу
        MATRIX res(columns_, rows_);
        for (int i = 0; i < rows_; ++i)
            for (int j = 0; j < columns_; ++j)
                res(i, j) = pow(-1, i + j) * Minordet(i, j);
        return res;
    }
    double ad(int row, int col) const
    { // алг доп для эл-та по его индексу
        return pow(-1, row + col) * Minordet(row, col);
    }
    MATRIX adMatrix() const
    { // матрица алг доп
        MATRIX res(rows_, columns_);
        for (int i = 0; i < rows_; ++i)
        {
            for (int j = 0; j < columns_; ++j)
            {
                res(i, j) = ad(i, j);
            }
        }
        return res;
    }
    MATRIX inv() const
    { // обратная матрица
        if (!isSqr())
            return MATRIX(0, 0);
        double deter = det();
        if (deter == 0)
            return MATRIX(0, 0);
        return adMatrix().transp() / deter;
    }

    friend MATRIX Solve(const MATRIX &A, const MATRIX &B);
    friend std::ostream &operator<<(std::ostream &out, const MATRIX &m);
};

MATRIX Solve(const MATRIX &A, const MATRIX &B)
{ // решение матричного ур-ния
    MATRIX InvA = A.inv();
    MATRIX res(InvA.rows_, B.columns_);
    for (int i = 0; i < InvA.rows_; ++i)
    {
        for (int j = 0; j < B.columns_; ++j)
        {
            res(i, j) = 0;
            for (int k = 0; k < InvA.columns_; ++k)
                res(i, j) += InvA(i, k) * B(k, j);
        }
    }
    return res;
}

std::ostream &operator<<(std::ostream &out, const MATRIX &m)
{
    for (const auto &row : m.data)
    {
        for (const auto &elem : row)
        {
            out << elem << "\t";
        }
        out << std::endl;
    }
    return out;
}

class Vector : public MATRIX
{
public:
    Vector() : MATRIX(1, 0) {}
    Vector(int size) : MATRIX(size, 1) {}
    Vector(const VECTOR &vect) : MATRIX(vect.size(), 1)
    {
        for (int i = 0; i < getRows(); i++)
        {
            data[i][0] = vect[i];
        }
    }
    Vector(const MATRIX &matr) : MATRIX(matr)
    {
        if (matr.getColumns() != 1)
            std::cout << "MATRIX must have one column." << std::endl;
    }
    ~Vector() {}

    double &operator()(int i) { return data[i][0]; }
    double operator()(int i) const { return data[i][0]; }
    Vector operator+(const Vector &matr) const { return static_cast<Vector>(MATRIX::operator+(matr)); }
    Vector operator-(const Vector &matr) const { return static_cast<Vector>(MATRIX::operator-(matr)); }
    Vector operator*(double digit) const { return static_cast<Vector>(MATRIX::operator*(digit)); }
    Vector operator/(double digit) const { return static_cast<Vector>(MATRIX::operator/(digit)); }
    bool operator==(const Vector &matr) const { return MATRIX::operator==(matr); }
    bool operator!=(const Vector &matr) const { return MATRIX::operator!=(matr); }
};

int main()
{
    std::vector<VECTOR> matr1 = {{3, 7}, {2, 8}};
    MATRIX A(matr1);
    std::cout << "MATRIX A:\n"
              << A << std::endl;
    std::cout << "det A:\t" << A.det() << std::endl
              << std::endl;
    std::cout << "transp matrix A:\n"
              << A.transp() << std::endl;
    std::cout << "Minor matrix A:\n"
              << A.MMatrix(0, 0) << std::endl;
    std::cout << "inv matrix A:\n"
              << A.inv() << std::endl;

    std::vector<VECTOR> matr2 = {{4, 8}, {6, 2}};
    MATRIX B(matr2);
    std::cout << "MATRIX B: " << std::endl;
    std::cout << B << std::endl;

    MATRIX X = Solve(A, B);
    std::cout << "Solving the matrix eguation AX=B: " << std::endl;
    std::cout << X << std::endl;

    VECTOR v = {1, 2, 3, 4, 5};
    Vector vect(v);
    std::cout << "Vector vect:\n"
              << vect << std::endl;
}