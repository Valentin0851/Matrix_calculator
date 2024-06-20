#include <iostream>
#include <vector>
#include <cmath>
#include <initializer_list>

template <class T>
using VECTOR = std::vector<T>;

template <class T>
using MASSIVE = std::vector<std::vector<T>>;

template <class T>
class MATRIX
{
    // Матрица организована как вектор строк(row),
    // т.е. сначала идет строка, потом обращение к колонке
protected:
    unsigned int rows_, columns_;
    MASSIVE<T> data;

public:
    MATRIX() = default;

    MATRIX(int rows, int columns)
    {
        this->rows_ = rows;
        this->columns_ = columns;
        this->data.resize(rows);
        for (auto &row : data)
        {
            row.resize(columns);
        }
    };
    MATRIX(const std::vector<std::vector<T>> matr)
    {
        data = matr;
        this->rows_ = matr.size();
        this->columns_ = matr[0].size();
    };

    MATRIX(const MATRIX &m)
    {
        this->rows_ = m.rows_;
        this->columns_ = m.columns_;
        for (int i = 0; i < m.rows_; i++)
        {
            for (int j = 0; j < m.columns_; j++)
            {
                this->data[i][j] = m.data[i][j];
            }
        }
    };

    ~MATRIX() = default;

    int getRowCount()
    {
        return this->rows_;
    };
    int rowLen() const
    {
        return this->columns_;
    };

    int getColCount()
    {
        return this->columns_;
    };
    int colLen() const
    {
        return this->rows_;
    };

    std::pair<int, int> dim() const
    {
        return {this->rows_, this->columns_};
    };

    void resize(int newRow, int newCol, bool &isSafeResize = true)
    {
        isSafeResize = (newRow > 0) && (newCol > 0);
        if (isSafeResize)
        {
            this->rows_ = newRow;
            this->columns_ = newCol;
            this->data.resize(newRow);
            for (auto &row : data)
            {
                row.resize(newCol);
            }
        }
    };

    VECTOR<T> toString() const
    {
        VECTOR<T> ans;
        for (int i = 0; i < this->rows_; i++)
        {
            for (int j = 0; j < this->columns_; j++)
            {
                ans.push_back(this->data[i][j]);
            }
        }
        return ans;
    };

    void print()
    {
        for (int i = 0; i < rows_; i++)
        {
            for (int j = 0; j < columns_; j++)
            {
                std::cout << this->data[i][j] << ' ';
            }
            std::cout << std::endl;
        }
    };

    void setItem(int row, int columns, T value)
    {
        data[row][columns] = value;
    };

    void setRow(int rNum, VECTOR<T> row)
    {
        for (int i = 0; i < columns_; i++)
        {
            data[rNum][i] = row[i];
        }
    };

    void setCol(int cNum, VECTOR<T> col)
    {
        for (int i = 0; i < rows_; i++)
        {
            data[i][cNum] = row[i];
        }
    };

    T it(int row, int columns) const
    {
        return data[row][columns];
    };

    VECTOR<T> row(int indx)
    {
        VECTOR<T> ans(columns_);
        for (int i = 0; i < columns_; i++)
        {
            ans[i] = data[indx][i];
        }
        return ans;
    };

    VECTOR<T> col(int indx)
    {
        VECTOR<T> ans(rows_);
        for (int i = 0; i < rows_; i++)
        {
            ans[i] = data[i][indx];
        }
        return ans;
    };

    const MASSIVE<T> &matr() const { return data; } // для константных объектов

    T det()
    {
        if (!this->isSqr())
        {
            return 0;
        }
        if (this->rows_ == 1)
        {
            return this->data[0][0];
        }
        if (this->rows_ == 2)
        {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }
        T ans = 0;

        for (int i = 0; i < this->columns_; i++)
        {
            ans += (i % 2 == 0 ? 1 : -1) * data[0][i] * Minor(0, i);
        }
        return ans;
    }; // Вычисление определителя

    MATRIX TRANS() // Транспонированная матрица
    {
        MATRIX ans(columns_, rows_);
        for (int i = 0; i < columns_; i++)
        {
            for (int j = 0; j < rows_; j++)
            {
                ans.data[i][j] = this->data[j][i];
            }
        }
        return ans;
    };

    MATRIX M(int row, int col) // Матрицу минора для элемента по его индексу
    {
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
    };

    T Minor(int row, int col)
    { // Возвращает определитель М. минора для элемента
        return this->M(row, col).det();
    };

    MATRIX MMatr()
    {
        MATRIX ans(this->rows_, this->columns_);
        for (int i = 0; i < this->rows_; i++)
        {
            for (int j = 0; j < this->columns_; j++)
            {
                ans.data[i][j] = this->Minor(i, j);
            }
        }
        return ans;
    }; // Матрица минора для элементов

    T ad(int row, int col)
    {
        return std::pow(-1, row + col) * Minor(row, col);
    }; // Алгебраическое дополнение
    MATRIX adMatr()
    {
        MATRIX ans(this->rows_, this->columns_);
        for (int i = 0; i < this->rows_; i++)
        {
            for (int j = 0; j < this->columns_; j++)
            {
                ans.data[i][j] = this->ad(i, j);
            }
        }
        return ans;
    }; // Матрица алгебраических дополнений

    bool isE()
    {
        bool ans = true;
        for (int i = 0; i < rows_; i++)
        {
            for (int j = 0; j < columns_; j++)
            {
                ans = ans && ((i == j) && (this->data[i][i] == 1) || (i != j) && (this->data[i][j] == 0));
            }
        }
        return ans;
    }; // Единичная ли матрица?
    bool isSqr()
    {
        return this->rows_ == this->columns_;
    }; // Квадратная ли матрица?

    

public:
    T &operator()(int row, int col)
    {
        return this->data[row][col];
    };
    const T &operator()(int row, int col) const
    {
        return this->data[row][col];
    }; // для константных объектов

    MATRIX &operator=(const MATRIX &m)
    {
        for (int i = 0; i < rows_; i++)
        {
            for (int j = 0; j < columns_; j++)
            {
                this->data[i][j] = m.data[i][j];
            }
        }
        return *this;
    };

    MATRIX operator+(const MATRIX &matr)
    {
        for (int i = 0; i < rows_; i++)
        {
            for (int j = 0; j < columns_; j++)
            {
                this->data[i][j] += matr[i][j];
            }
        }
        return *this;
    };
    MATRIX operator-(const MATRIX &matr)
    {
        for (int i = 0; i < rows_; i++)
        {
            for (int j = 0; j < columns_; j++)
            {
                this->data[i][j] -= matr[i][j];
            }
        }
        return *this;
    };
    MATRIX operator*(const double &l)
    {
        for (int i = 0; i < rows_; i++)
        {
            for (int j = 0; j < columns_; j++)
            {
                this->data[i][j] *= l;
            }
        }
        return *this;
    };
    MATRIX operator/(const double &l)
    {
        MATRIX res(rows_, columns_);
        for(int i = 0; i < rows_; ++i)
            for(int j = 0; j < columns_; ++j)
                res.data[i][j] = this->data[i][j]/l;
        return res;
    };

    friend MATRIX operator*(const MATRIX &Matr1, const MATRIX &Matr2)
    {
        MATRIX ans(Matr1.rows_, Matr2.columns_);
        for (int i = 0; i < Matr1.rows_; i++)
        {
            for (int j = 0; j < Matr2.columns_; j++)
            {
                ans.data[i][j] = 0;

                for (int k = 0; k < Matr1.columns_; k++)
                {
                    ans.data[i][j] += Matr1.data[i][k] * Matr2.data[k][j];
                }
            }
        }
        return ans;
    };
    friend MATRIX operator*(const double num, const MATRIX &Matr)
    {
        for (int i = 0; i < Matr.getRowCount(); i++)
        {
            for (int j = 0; j < Matr.getColCount(); j++)
            {
                Matr.data[i][j] *= num;
            }
        }
        return Matr;
    };

    friend bool operator==(const MATRIX &m1, const MATRIX &m2)
    {
        bool ans = (m1.dim().first == m2.dim().first) && (m1.dim().second == m2.dim().second);
        if (ans)
        {
            for (int i = 0; i < m1.getRowCount(); i++)
            {
                for (int j = 0; j < m1.getColCount(); j++)
                {
                    ans = ans && m1.data[i][j] == m2.data[i][j];
                }
            }
        }
        return ans;
    };

    friend bool operator!=(const MATRIX &m1, const MATRIX &m2)
    {
        bool ans = (m1.dim().first == m2.dim().first) && (m1.dim().second == m2.dim().second);
        if (ans)
        {
            for (int i = 0; i < m1.getRowCount(); i++)
            {
                for (int j = 0; j < m1.getColCount(); j++)
                {
                    ans = ans && m1.data[i][j] == m2.data[i][j];
                }
            }
        }
        if (ans)
        {
            return false;
        }
        else
        {
            return true;
        }
    };

    friend std::istream &operator>>(std::istream &in, MATRIX &m)
    {
        for (const auto &row : m.data)
        {
            for (const auto &elem : row)
            {
                in >> elem << "\t";
            }
        }
        return in;
    };

    friend std::ostream &operator<<(std::ostream &out, const MATRIX &m)
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
    };
    
    
    MATRIX inv()
    {
        if (!isSqr())
        {
            return MATRIX(0, 0);
        }
        T determinate = this->det();
        if (determinate == 0)
        {
            return MATRIX(0, 0);
        }
        return this->adMatr().TRANS() / determinate;
    }; // Обратная матрица

    static MATRIX solve(const MATRIX &A, const MATRIX &B)
    {
        MATRIX invA = A.inv();
        MATRIX ans(invA.getRowCount(), B.getColCount());
        for (int i = 0; i < invA.getRowCount(); i++)
        {
            for (int j = 0; j < B.getColCount(); j++)
            {
                ans.data[i][j] = 0;
                for (int k = 0; k < invA.getColCount(); k++)
                {
                    ans.data[i][j] += invA.data[i][k] * B.data[k][j];
                }
            }
        }
        return ans;
    }; // Решение матричного уравнения A*X=B
};

template <class T>
class Vector : public MATRIX<T>
{
public:
    Vector() : MATRIX<T>(1, 0) {}
    Vector(int size) : MATRIX<T>(size, 1) {}
    Vector(const std::vector<T> &vect) : MATRIX<T>(vect.size(), 1)
    {
        for (int i = 0; i < Vector<T>::getRowCount(); i++)
        {
            this->data[i][0] = vect[i];
        }
    }
    Vector(const MATRIX<T> &other) : MATRIX<T>(other)
    {
        if (other.getColumnCount() != 1)
            std::cout << "Matrix must have one column." << std::endl;
    }
    ~Vector() {}

    double &operator()(int i) { return this->data[i][0]; }
    double operator()(int i) const { return this->data[i][0]; }
    Vector operator+(const Vector &other) const { return static_cast<Vector>(MATRIX<T>::operator+(other)); }
    Vector operator-(const Vector &other) const { return static_cast<Vector>(MATRIX<T>::operator-(other)); }
    Vector operator*(double scalar) const { return static_cast<Vector>(MATRIX<T>::operator*(scalar)); }
    Vector operator/(double scalar) const { return static_cast<Vector>(MATRIX<T>::operator/(scalar)); }
    bool operator==(const Vector &other) const { return MATRIX<T>::operator==(other); }
    bool operator!=(const Vector &other) const { return MATRIX<T>::operator!=(other); }
};

int main()
{
    std::vector<std::vector<double>> a = {{1, 3, 5}, {2, 4, 6}, {7, 8, 10}};
    MATRIX<double> A(a);

    std::cout << "Start Matrix: \n";
    A.print();
    std::cout << "\n";

    std::cout << "Transopned matrix: \n";
    A.TRANS().print();
    std::cout << "\n";

    std::cout << "Determinant A:\t" << A.det() << std::endl
              << std::endl;

    std::cout << "MMatrix: \n";
    A.MMatr().print();

    std::cout << "alg add: \n";
    A.adMatr().print();
    
    std::cout<< "inverse matrix: \n";
    A.inv().print();
    std::cout<<std::endl;
    MATRIX<double> B = A.inv();

    std::cout<<"proverka\n";
    MATRIX<double> C = A*B;
    C.print();
    return 0;
}