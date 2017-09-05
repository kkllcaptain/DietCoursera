/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package diet;

/**
 *
 * @author sguo46
 */
import java.io.*;
import java.util.*;

public class Diet {

    BufferedReader br;
    PrintWriter out;
    StringTokenizer st;
    boolean eof;

    int solveDietProblem(int n, int m, double A[][], double[] b, double[] c, double[] x) {
        Arrays.fill(x, 1);
        // Write your code here - simplex algorithm
        
        // manage the status
        int status = 0; 
        
        // permutate all m out of n euqations to construct a system of euqalities to solve for a vetex
        List<List> Sets = new ArrayList<>();
        
        Sets = chooseSubSet(n,m);
        if((int)Sets.get(0).get(0) == -1) status = 1;
        status = processStatus(status, A, b);
        if (status==-1 || status ==1) return status;

        
        // now we have a subsets, lets print out all sets
        for(int setIter = 0; setIter < Sets.size(); setIter ++){
            List subset = Sets.get(setIter);
            
            System.out.println(Arrays.toString(subset.toArray()));

            
            
            
            
        }
        
        
        
        
        // generate valid subsets
        double[][] A_E = A;
        double[] b_E = b;
        
        // check if the selected constraints are parallel, if so, skip, and output an warning that there might not be an answer
        
        // solve the system of euqalities by Gaussian elimination       
        Equation newEq = new Equation(A_E, b_E);
        
        // the obtained is a vertex
        
        // relax one constraint to obtain an edge (that connect to the vertex)
        
        // return -1 if infinite
        
        // find the other vertex on the edge
        
        

        return 0;
    }
    
    static int processStatus(int status, double A[][], double[] b){
        //if(status == 1){ // not possible for bounded solution, still possible to have no solution
            for(int row =0; row<A.length; row++){
                for(int otherRow = 0; otherRow<A.length; otherRow ++){
                    if(row != otherRow){
                        double ratio = A[row][0]/A[otherRow][0];
                        boolean parallel = true;
                        for(int column = 1; column <A[0].length; column++){
                            if(A[row][column]/A[otherRow][column] != ratio){
                                parallel = false;
                            }
                        }
                        if(parallel&&(b[row] != b[otherRow])){
                            status = -1;
                            return status;
                        } else return status;
                    }
                }
            }
        //}
        
        return status;
        
    }
    
    static List<List> chooseSubSet(int n, int m){
        // choose n euqation out of m equations
        
        List<List> Sets = new ArrayList<>();
        // if there are more variables than constraints
        if (m>n){
            List<Integer> subset = new ArrayList<>();
            subset.add(-1);
            Sets.add(subset);
            return Sets;
        } else if(m == n){
            List<Integer> subset = new ArrayList<>();
            for (int iter = 0; iter<m; iter ++){
                subset.add(1);
            }
            Sets.add(subset);
            return Sets;
        } else {

            for (int i = 0; i < (1 << n); i++) {
                ArrayList<Integer> subset = new ArrayList<>();
                for (int j = 0; j < n; j++) {
                    if (((i >> j) & 1) == 1) {
                        subset.add(1);
                    } else subset.add(0);
                }
                
                int sum = subset.stream().mapToInt(Integer::intValue).sum();
                if (sum == m) Sets.add(subset);                
                
            }

            return Sets;
        }
        
        

        
        
    }

    static double[] SolveEquation(Equation equation) {
        double a[][] = equation.a;
        double b[] = equation.b;
        int size = a.length;

        boolean[] used_columns = new boolean[size];
        boolean[] used_rows = new boolean[size];
        for (int step = 0; step < size; ++step) {
            Position pivot_element = SelectPivotElement(a, used_rows, used_columns);
            SwapLines(a, b, used_rows, pivot_element);
            ProcessPivotElement(a, b, pivot_element);
            MarkPivotElementUsed(pivot_element, used_rows, used_columns);
        }

        return b;
    }

    static Position SelectPivotElement(double a[][], boolean used_rows[], boolean used_columns[]) {
        // This algorithm selects the first free element.
        // You'll need to improve it to pass the problem.
        Position pivot_element = null;
        pivot_element.column =0;
        pivot_element.row =0;

        while (used_rows[pivot_element.row]) {
            ++pivot_element.row;
        }
        while (used_columns[pivot_element.column] || a[pivot_element.row][pivot_element.column] == 0) {

            ++pivot_element.column;
            if (pivot_element.column >= a[0].length) {
                pivot_element.row = -1;
                pivot_element.column = -1;
                break;
            }
        }

        return pivot_element;
    }

    static void SwapLines(double a[][], double b[], boolean used_rows[], Position pivot_element) {
        int size = a.length;

        if (pivot_element.row == -1) {
            throw new IllegalArgumentException("\n\n ====== \nNo valid Pivot for this operation\n ====== \n");
        }

        for (int column = 0; column < size; ++column) {
            double tmpa = a[pivot_element.column][column];
            a[pivot_element.column][column] = a[pivot_element.row][column];
            a[pivot_element.row][column] = tmpa;
        }

        double tmpb = b[pivot_element.column];
        b[pivot_element.column] = b[pivot_element.row];
        b[pivot_element.row] = tmpb;

        boolean tmpu = used_rows[pivot_element.column];
        used_rows[pivot_element.column] = used_rows[pivot_element.row];
        used_rows[pivot_element.row] = tmpu;

        pivot_element.row = pivot_element.column;
    }

    static void ProcessPivotElement(double a[][], double b[], Position pivot_element) {
        // Write your code here
        int pivot_row = pivot_element.row;
        int pivot_column = pivot_element.column;

        if (pivot_row != -1) {
            double pivot_value = a[pivot_row][pivot_column];

//        System.out.println("=========== Iteration ===============");
//        for (int r=0; r<a.length; r++){
//            for (int c=0; c<a[0].length; c++){
//                System.out.printf("%f ", a[r][c]);
//            }
//            System.out.print(b[r]);
//            System.out.println();
//        }
//        
            // process the pivot row
            for (int c = 0; c < a[0].length; c++) {
                a[pivot_row][c] /= (double) pivot_value;
            }
            b[pivot_row] /= (double) pivot_value;

            //System.out.printf("the pivot was chosen to be the element of (%s, %s) with the value of %f \n", pivot_row, pivot_column, a[pivot_row][pivot_column]);
            // process other row
            for (int r = 0; r < a.length; r++) {
                if (r != pivot_row) {
                    double ratio = a[r][pivot_column] / 1.0;

                    for (int c = 0; c < a[0].length; c++) {
                        a[r][c] -= ratio * a[pivot_row][c];
                    }

                    b[r] -= ratio * b[pivot_row];
                }
            }

//        System.out.println("=========== after Process ===============");
//        for (int r=0; r<a.length; r++){
//            for (int c=0; c<a[0].length; c++){
//                System.out.printf("%f ", a[r][c]);
//            }
//            System.out.print(b[r]);
//            System.out.println();
//        }
        } else {
            throw new IllegalArgumentException("No valid Pivot for this operation");
        }
    }

    static void MarkPivotElementUsed(Position pivot_element, boolean used_rows[], boolean used_columns[]) {
        used_rows[pivot_element.row] = true;
        used_columns[pivot_element.column] = true;
    }

    class Equation {

        Equation(double a[][], double b[]) {
            this.a = a;
            this.b = b;
        }

        double a[][];
        double b[];
    }

    class Position {

        Position(int column, int row) {
            this.column = column;
            this.row = row;
        }

        int column;
        int row;
    }

    void solve() throws IOException {
        int n = nextInt(); // number of restrictions/constraints
        int m = nextInt(); // number of variables

        // standard linear in equalities "Ax <= b"
        // A is nXm matrix, x is a vector of mX1 and b is vector of nX1
        double[][] A = new double[n][m];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                A[i][j] = nextInt();
            }
        }
        double[] b = new double[n];
        for (int i = 0; i < n; i++) {
            b[i] = nextInt();
        }

        // weights of objective function of each different items, c is a vector of size m
        double[] c = new double[m];
        for (int i = 0; i < m; i++) {
            c[i] = nextInt();
        }

        // solve for x
        double[] ansx = new double[m];
        int anst = solveDietProblem(n, m, A, b, c, ansx);
        if (anst == -1) {
            out.printf("No solution\n");
            return;
        }
        if (anst == 0) {
            out.printf("Bounded solution\n");
            for (int i = 0; i < m; i++) {
                out.printf("%.18f%c", ansx[i], i + 1 == m ? '\n' : ' ');
            }
            return;
        }
        if (anst == 1) {
            out.printf("Infinity\n");
            return;
        }
    }

    Diet() throws IOException {
        br = new BufferedReader(new InputStreamReader(System.in));
        out = new PrintWriter(System.out);
        solve();
        out.close();
    }

    public static void main(String[] args) throws IOException {
        new Diet();
    }

    String nextToken() {
        while (st == null || !st.hasMoreTokens()) {
            try {
                st = new StringTokenizer(br.readLine());
            } catch (Exception e) {
                eof = true;
                return null;
            }
        }
        return st.nextToken();
    }

    int nextInt() throws IOException {
        return Integer.parseInt(nextToken());
    }

}
