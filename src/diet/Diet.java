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
        // m restrictions
        // n variables
        
        // ====================================================
        // [Preprocessing] remove all parallel euqations
        // if incorrect euqations were found, return as no solution (-1)
        // if more variables than restrictions, return as infinit (+1)
        // ====================================================
        
        List<ArrayList> AList = new ArrayList<>();
        List<Double> bList = new ArrayList<>();
        for(int i =0; i<A.length; i++){
            List RowList = new ArrayList<Double>();
            for(int j=0; j<A[0].length; j++){
                RowList.add(A[i][j]);
            }
            AList.add(new ArrayList(RowList));
            bList.add(b[i]);
        }
        
        int newSize = preProcessEquations(AList, bList);
        if (newSize == -1) return -1;
        else n -= newSize;
        //System.out.println("the number of equations that were removed: " + newSize);
        if (n<m) return 1;


        
        // ================= End of PreProcessing =========================
        
        //  
        // permutate all m out of n euqations to construct a system of euqalities to solve for a vetex
        List<List> Sets = new ArrayList<>();
        Sets = chooseSubSet(n,m);
                
        // now we have a subsets, lets print out all sets
        //for(int setIter = 0; setIter < Sets.size(); setIter ++){
        List subset = Sets.get(0);// may not need iterate all sets. 
  
        Equation newEq = constructNewEquations(AList, bList, subset);
     
        
            
        System.out.println(Arrays.toString(subset.toArray()));


        // solve the system of euqalities by Gaussian elimination       
        x = SolveEquation(newEq);
        
        // the obtained is a vertex
        
        // relax one constraint to obtain an edge (that connect to the vertex)
        
        // return -1 if infinite
        
        // find the other vertex on the edge
        
        

        return 0;
    }
    
    int preProcessEquations(List<ArrayList> A, List<Double> b){
        int removed = 0; // no restriction was removed (yet)
        for(int row =0; row<A.size(); row++){
                for(int otherRow = 0; otherRow<A.size(); otherRow++){
                    if(row != otherRow){
                        double ratio = 0.0;
                        int cc = 0;
                        
                        while (ratio == 0.0 && cc < A.size()){
                            if( (double) A.get(row).get(cc) != 0) {
                                ratio = ((double)A.get(row).get(cc))/((double) A.get(otherRow).get(cc));
                            } else {
                                cc++;
                            }
                            }
                        boolean Aparallel = true;
                        for(int column = 0; column <A.get(0).size(); column++){
                            double currentRatio = 0.0;
                            if (Objects.equals((double)A.get(row).get(column), (double) A.get(otherRow).get(column))&&(double)A.get(row).get(column) == 0.0) {
                                currentRatio = ratio;
                            } else {
                                currentRatio = ((double)A.get(row).get(column))/((double) A.get(otherRow).get(column));
                            }
                            //System.out.println("the current raio is " + currentRatio + " and the ratio is " + ratio);
                            if(!Objects.equals(currentRatio, ratio)){
                                Aparallel = false;
                                // not parallel
                            }
                        }
                        boolean bParallel = (Objects.equals(b.get(row), b.get(otherRow)));
                        //System.out.println("Row "+row+" and Row "+otherRow+" is "+Aparallel+" and "+bParallel);
                        if(Aparallel&&bParallel){
                            // truely parallel
                            A.remove(otherRow);
                            b.remove(otherRow);
                            removed++;
                        } else if (Aparallel&&(!bParallel)){
                            // incorrect fomulation
                            A.clear();
                            b.clear();
                            return -1;
                        }
                    }
                }
            }
        return removed;
    }
    
    Equation constructNewEquations(List<ArrayList> A, List<Double> b, List subset){
        List<ArrayList> AE = A;
        List<Double> bE = b;
        
        for(int index =0; index<A.size(); index++){
            if((int)subset.get(index) == 0){
                AE.remove(index);
                bE.remove(index);
            }
        }

        double[][] Ap = new double[AE.size()][AE.get(0).size()];
        double[] bp = new double[AE.size()];
        
        for (int row = 0; row < AE.size(); row++){
            bp[row] = (double) bE.get(row);
            for (int column = 0; column < AE.get(0).size(); column ++){
                Ap[row][column] = (double) AE.get(row).get(column);
            }
            
        }
        
        Equation myEq = new Equation(Ap,bp);
        
        return myEq;
        
    }

    List<List> chooseSubSet(int n, int m){
        // choose n euqation out of m equations
        
        List<List> Sets = new ArrayList<>();
        // if there are more variables than constraints
        if(m == n){
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

    double[] SolveEquation(Equation equation) {
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

    Position SelectPivotElement(double a[][], boolean used_rows[], boolean used_columns[]) {
   
        Position pivot_element = new Position(0,0);

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

    void SwapLines(double a[][], double b[], boolean used_rows[], Position pivot_element) {
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

    void ProcessPivotElement(double a[][], double b[], Position pivot_element) {
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

    void MarkPivotElementUsed(Position pivot_element, boolean used_rows[], boolean used_columns[]) {
        used_rows[pivot_element.row] = true;
        used_columns[pivot_element.column] = true;
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