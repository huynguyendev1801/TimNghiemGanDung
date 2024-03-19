#include <iostream>
#include <cmath>
#include <fstream>
#include <conio.h>   // For arrow key input
#include <windows.h> // For system("cls") and Sleep()
using namespace std;

// Tinh chuan cua ma tran
double matrix_norm(double A[][100], int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        double row_sum = 0;
        for (int j = 0; j < n; j++)
        {
            norm += pow(A[i][j], 2);
        }
    }
    norm = sqrt(norm);
    return norm;
}

double vector_norm(double x[], int n)
{
    double norm = 0;
    for (int i = 0; i < n; i++)
    {
        norm += x[i] * x[i];
    }
    return sqrt(norm);
}

// Phuong phap lap don
void jacobi_iter(double A[][100], double B[], double x[], int n, double tol, int max_iter)
{
    ofstream output_file("output3.txt", ios::app); // Mở file ở chế độ append
    double D[100][100], L[100][100], U[100][100];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                D[i][j] = A[i][j];
                L[i][j] = 0;
                U[i][j] = 0;
            }
            else if (i < j)
            {
                L[i][j] = A[i][j];
                D[i][j] = 0;
                U[i][j] = 0;
            }
            else
            {
                L[i][j] = 0;
                D[i][j] = 0;
                U[i][j] = A[i][j];
            }
        }
    }
    // Thuc hien phuong phap lap don
    int iter = 0;
    double x_new[n];
    double diff = tol + 1;
    while (diff > tol && iter < max_iter)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (1 / A[i][i]) * (B[i] - sum);
        }

        diff = 0;
        for (int i = 0; i < n; i++)
        {
            diff += pow(x_new[i] - x[i], 2);
        }
        diff = sqrt(diff);

        for (int i = 0; i < n; i++)
        {
            x[i] = x_new[i];
        }

        iter++;
    }
    // In ket qua
    if (iter >= max_iter)
    {
        output_file << "Phuong phap lap don khong hoi tu sau " << max_iter << " vong lap." << endl;
    }
    else
    {
        output_file << "Phuong phap lap don da hoi tu sau " << iter << " vong lap." << endl;
        output_file << "Nghiem cua he phuong trinh la: ";
        for (int i = 0; i < n; i++)
        {
            output_file << x[i] << " ";
        }
        output_file << endl;
    }
    output_file.close();
}

// Phuong phap lap Seidel
void gauss_seidel_iter(double A[][100], double B[], double x[], int n, double tol, int max_iter)
{
    ofstream output_file("output3.txt", ios::app); // Mở file ở chế độ append
    double D[100][100], L[100][100], U[100][100];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                D[i][j] = A[i][j];
                L[i][j] = 0;
                U[i][j] = 0;
            }
            else if (i < j)
            {
                L[i][j] = A[i][j];
                D[i][j] = 0;
                U[i][j] = 0;
            }
            else
            {
                L[i][j] = 0;
                D[i][j] = 0;
                U[i][j] = A[i][j];
            }
        }
    }
    // Thuc hien phuong phap lap Seidel
    int iter = 0;
    double x_new[n];
    double diff = tol + 1;
    while (diff > tol && iter < max_iter)
    {
        // Tinh vecto x moi
        for (int i = 0; i < n; i++)
        {
            x_new[i] = 0;
            for (int j = 0; j < n; j++)
            {
                if (j < i)
                {
                    x_new[i] += -L[i][j] * x_new[j];
                }
                else if (j > i)
                {
                    x_new[i] += -U[i][j] * x[j];
                }
            }
            x_new[i] /= D[i][i];
            x_new[i] += B[i] / D[i][i];
        }

        // Tinh 2 so giua 2 vecto x lien tiep
        diff = 0;
        for (int i = 0; i < n; i++)
        {
            diff += pow(x_new[i] - x[i], 2);
        }
        diff = sqrt(diff);

        // Cap nhat vector x
        for (int i = 0; i < n; i++)
        {
            x[i] = x_new[i];
        }

        iter++;
    }

    // In ket qua
    if (iter >= max_iter)
    {
        output_file << "Phuong phap lap Seidel khong hoi tu sau " << max_iter << " vong lap." << endl;
    }
    else
    {
        output_file << "Phuong phap lap Seidel da hoi tu sau " << iter << " vong lap." << endl;
        output_file << "Nghiem cua he phuong trinh la: ";
        for (int i = 0; i < n; i++)
        {
            output_file << x[i] << " ";
        }
        output_file << endl;
    }
    output_file.close();
}

double get_spectral_radius_jacobi(double A[100][100], int N)
{
    ofstream output_file("output2.txt", ios::app); // Mở file ở chế độ append
    double D[100][100], R[100][100];
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                D[i][j] = A[i][j];
                R[i][j] = 0;
            }
            else
            {
                D[i][j] = 0;
                R[i][j] = A[i][j];
            }
        }
    }
    double B[N][N];
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            B[i][j] = -D[i][j] / A[i][i];
        }
    }
    double eigenvalues[N];
    int max_iter = 100;      // so lan lap toi da
    double tol = 1e-6;       // sai so tu do cho phep
    double x[N] = {1, 1, 1}; // vecto khoi tao
    double x_old[N];
    for (i = 0; i < max_iter; i++)
    {
        for (j = 0; j < N; j++)
        {
            x_old[j] = x[j];
            double sum = 0;
            int k;
            for (k = 0; k < N; k++)
            {
                sum += B[j][k] * x[k];
            }
            x[j] = sum;
        }
        double err = 0;
        for (j = 0; j < N; j++)
        {
            err += fabs(x[j] - x_old[j]);
        }
        if (err < tol)
        {
            break;
        }
    }
    if (i == max_iter)
    {
        output_file << "Khong tim thay gia tri rieng lon nhat.\n";
        return -1;
    }
    for (i = 0; i < N; i++)
    {
        eigenvalues[i] = x[i] / x_old[i];
    }
    double max_eigenvalue = eigenvalues[0];
    for (i = 1; i < N; i++)
    {
        if (fabs(eigenvalues[i]) > fabs(max_eigenvalue))
        {
            max_eigenvalue = eigenvalues[i];
        }
    }
    output_file.close();
    return fabs(max_eigenvalue);
}

double get_spectral_radius_seidel(double A[100][100], int N)
{
    ofstream output_file("output2.txt", ios::app); // Mở file ở chế độ append
    double L[100][100], D[100][100], U[100][100];
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i < j)
            {
                U[i][j] = A[i][j];
                L[i][j] = 0;
                D[i][j] = 0;
            }
            else if (i == j)
            {
                D[i][j] = A[i][j];
                U[i][j] = 0;
                L[i][j] = 0;
            }
            else
            {
                L[i][j] = A[i][j];
                D[i][j] = 0;
                U[i][j] = 0;
            }
        }
    }
    double B[N][N];
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i == j)
            {
                B[i][j] = -D[i][j] / (L[i][j] + D[i][j]);
            }
            else
            {
                B[i][j] = -(L[i][j] / (L[i][i] + D[i][i])) - (U[i][j] / (D[j][j] + U[j][j]));
            }
        }
    }
    double eigenvalues[N];
    int max_iter = 100;      // so lan lap toi da
    double tol = 1e-6;       // sai so tu do cho phep
    double x[N] = {1, 1, 1}; // vecto khoi tao
    double x_old[N];
    for (i = 0; i < max_iter; i++)
    {
        for (j = 0; j < N; j++)
        {
            x_old[j] = x[j];
            double sum = 0;
            int k;
            for (k = 0; k < N; k++)
            {
                sum += B[j][k] * x[k];
            }
            x[j] = sum;
        }
        double err = 0;
        for (j = 0; j < N; j++)
        {
            err += fabs(x[j] - x_old[j]);
        }
        if (err < tol)
        {
            break;
        }
    }
    if (i == max_iter)
    {
        output_file << "Khong tim thay gia tri rieng lon nhat.\n";
        return -1;
    }
    for (i = 0; i < N; i++)
    {
        eigenvalues[i] = x[i] / x_old[i];
    }
    double max_eigenvalue = eigenvalues[0];
    for (i = 1; i < N; i++)
    {
        if (fabs(eigenvalues[i]) > fabs(max_eigenvalue))
        {
            max_eigenvalue = eigenvalues[i];
        }
    }
    output_file.close();
    return fabs(max_eigenvalue);
}
// Phuong phap lap don
void jacobi_iterxk(double A[][100], double B[], double x[], int n, double tol, int max_iter)
{
    ofstream output_file("output4.txt", ios::app); // Mở file ở chế độ append
    double D[100][100], L[100][100], U[100][100];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                D[i][j] = A[i][j];
                L[i][j] = 0;
                U[i][j] = 0;
            }
            else if (i < j)
            {
                L[i][j] = A[i][j];
                D[i][j] = 0;
                U[i][j] = 0;
            }
            else
            {
                L[i][j] = 0;
                D[i][j] = 0;
                U[i][j] = A[i][j];
            }
        }
    }
    // Thuc hien phuong phap lap don
    int iter = 0;
    double x_new[100];
    double diff = tol + 1;
    while (diff > tol && iter < max_iter)
    {
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < n; j++)
            {
                if (j != i)
                {
                    sum += A[i][j] * x[j];
                }
            }
            x_new[i] = (1 / A[i][i]) * (B[i] - sum);
        }

        diff = 0;
        for (int i = 0; i < n; i++)
        {
            diff += fabs(x_new[i] - x[i]);
        }

        for (int i = 0; i < n; i++)
        {
            x[i] = x_new[i];
        }

        // Kiểm tra điều kiện hội tụ abs(X(k) - X(k-1)) <= e
        if (iter > 0 && diff <= tol)
        {
            break;
        }

        iter++;
    }

    // In ket qua
    if (iter >= max_iter)
    {
        output_file << "Phuong phap lap don khong hoi tu sau " << max_iter << " vong  voi X(k) thoa man: ||X(k) - X(k-1)|| <= e cho truoc." << endl;
    }
    else
    {
        output_file << "Phuong phap lap don da hoi tu sau " << iter << " vong lap voi X(k) thoa man: ||X(k) - X(k-1)|| <= e cho truoc." << endl;
        output_file << "Nghiem cua he phuong trinh la: ";
        for (int i = 0; i < n; i++)
        {
            output_file << x[i] << " ";
        }
        output_file << endl;
    }
    output_file.close();
}

int main()
{
    int n, max_iter;
    double tol;
    double A[100][100], B[100], x[100];
    bool is_diag_dominant = true;
    int selected_option = 1;

    while (true)
    {
        system("cls");

        if (1 == selected_option)
            cout << "->";
        else
            cout << "  ";
        cout << "Nhap vao A, b theo khuon dang cua ma tran" << endl;

        if (2 == selected_option)
            cout << "->";
        else
            cout << "  ";
        cout << "Kiem tra tinh cheo troi cua ma tran A" << endl;

        if (3 == selected_option)
            cout << "->";
        else
            cout << "  ";
        cout << "Tinh chuan cua ma tran A va kiem tra su hoi tu cua phuong phap lap don, lap Seidel doi voi he da cho" << endl;

        if (4 == selected_option)
            cout << "->";
        else
            cout << "  ";
        cout << "Tinh nghiem gan dung voi so lan lap k cho truoc va danh gia sai so, theo ca 2 cong thuc (k duoc nhap vao tu ban phim). Tinh nghiem gan dung voi so so e cho truoc, theo ca 2 cach ap dung cong thuc sai so (e duoc nhap tu ban phim)" << endl;

        if (5 == selected_option)
            cout << "->";
        else
            cout << "  ";
        cout << "Tinh nghiem gan dung X(k) thoa man: ||X(k) - X(k-1)|| <= e cho truoc" << endl;

        if (0 == selected_option)
            cout << "->";
        else
            cout << "  ";
        cout << "Thoat" << endl;
        // Arrow key navigation
        char key = getch();
        if (key == 13) // Enter key
        {
            // Perform action based on selected_option
             switch (selected_option)
        {
        case 1:
        {
            /* code */
            cout << "Nhap so chieu ma tran A: ";
            cin >> n;
            cout << "Nhap ma tran A " << n << "x" << n << ":" << endl;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    cin >> A[i][j];
                    // input_file >> A[i][j];
                }
            }
            cout << "Nhap vector B: " << endl;
            for (int i = 0; i < n; i++)
            {
                cin >> B[i];
                // input_file >> B[i]; chuyển thành menu điều hướng theo dạng nhân nút mũi ten lên xuống, 
            }
            // Nhap so lan lap toi da va sai so cho phep
            cout << "Nhap so lan lap toi da: ";
            cin >> max_iter;
            // input_file >> max_iter;
            cout << "Nhap sai so cho phep: ";
            cin >> tol;
            // input_file >> tol;
            // input_file.close();
        }

        break;
        case 2:
        {
            ofstream output_file("output.txt", ios::app); // Mở file ở chế độ append
                                                          // 2. Kiem tra tinh cheo troi cua ma tran A

            for (int i = 0; i < n; i++)
            {
                double diag = abs(A[i][i]);
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += abs(A[i][j]);
                    }
                }
                if (diag <= sum)
                {
                    is_diag_dominant = false;
                    break;
                }
            }

            if (is_diag_dominant)
                output_file << "Ma tran A la ma tran cheo troi\n";
            else
                output_file << "Khong the ap dung phuong phap lap don hoac phuong phap lap Seidel vi ma tran A khong cheo troi." << endl;
            output_file.close();
        }
        break;
        case 3:
        {                                                 // 2. Kiem tra tinh cheo troi cua ma tran A
            ofstream output_file("output2.txt", ios::app); // Mở file ở chế độ append
            for (int i = 0; i < n; i++)
            {
                double diag = abs(A[i][i]);
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += abs(A[i][j]);
                    }
                }
                if (diag <= sum)
                {
                    is_diag_dominant = false;
                    break;
                }
            }
            if (is_diag_dominant)
            {
                output_file << "Ma tran A la ma tran cheo troi\n";
                /* 3. Tinh chuan cua ma tran A va kiem tra su hoi tu cua phuong phap lap don,
                         lap Seidel doi voi he da cho */

                // 3.1: Tinh chuan cua ma tran A
                double norm_A = matrix_norm(A, n);
                output_file << "Chuan cua ma tran A: " << norm_A << endl;
                // 3.2: Kiem tra tinh hoi tu cua phuong phap lap don
                double spectral_radius_jacobi = get_spectral_radius_jacobi(A, n);
                if (spectral_radius_jacobi >= 1)
                    output_file << "Phuong phap lap don khong hoi tu" << endl;

                // 3.3: Kiem tra tinh hoi tu cua phuong phap lap Seidel
                double spectral_radius_seidel = get_spectral_radius_seidel(A, n);
                if (spectral_radius_seidel >= 1)
                    output_file << "Phuong phap lap Seidel khong hoi tu" << endl;
            }
            else
            {
                output_file << "Khong the ap dung phuong phap lap don hoac phuong phap lap Seidel vi ma tran A khong cheo troi." << endl;
            }
            output_file.close();
        }
        break;
        case 4:
        {
            // 2. Kiem tra tinh cheo troi cua ma tran A
            ofstream output_file("output3.txt", ios::app); // Mở file ở chế độ append
            for (int i = 0; i < n; i++)
            {
                double diag = abs(A[i][i]);
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += abs(A[i][j]);
                    }
                }
                if (diag <= sum)
                {
                    is_diag_dominant = false;
                    break;
                }
            }
            if (is_diag_dominant)
            {
                output_file << "Ma tran A la ma tran cheo troi\n";
                                      
                cout << "Nhap k:" << endl;
                int k; cin >>k;
                cout << "Nhap e:" << endl;
                double e;
                cin >> e;
                // 4,5: Tinh nghiem gan dung su dung phuong phap lap don
                for (int i = 0; i < n; i++)
                {
                    x[i] = 0;
                }
                jacobi_iter(A, B, x, n, e, k);

                // 4,5: Tinh nghiem gan dung su dung phuong phap lap Seidel
                for (int i = 0; i < n; i++)
                {
                    x[i] = 0;
                }
                gauss_seidel_iter(A, B, x, n, e, k);
            }
            else
            {
                output_file << "Khong the ap dung phuong phap lap don hoac phuong phap lap Seidel vi ma tran A khong cheo troi." << endl;
            }
            output_file.close();
        }
        break;
        case 5:
            // 6.
            for (int i = 0; i < n; i++)
            {
                x[i] = 0; // Khởi tạo nghiệm gần đúng x ban đầu
            }
            jacobi_iterxk(A, B, x, n, tol, max_iter);
            break;
        case 0:
            return 0;
        default:
            break;
        }
        }
        else if (key == 72) 
        {
            selected_option = (selected_option == 0) ? 5 : selected_option - 1;
        }
        else if (key == 80) 
        {
            selected_option = (selected_option == 5) ? 0 : selected_option + 1;
        }
        else
        {
            // Handle other keys
        }



    }

    return 0;
}