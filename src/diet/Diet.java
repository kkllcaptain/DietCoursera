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


public class Diet {

    BufferedReader br;
    PrintWriter out;
    StringTokenizer st;
    boolean eof;

    int solveDietProblem(int n, int m, double A[][], double[] b, double[] c, double[] x) {
        Arrays.fill(x, 1);
        // Write your code here
        
        // =============== pre-process ======================================
        // take matrix A and vector b, build list A and List b
        // add non-negtivity constraints
        // remove parallel constraints
        // detect incorrect (non-overlapping) conditions and return -1
        // modify Alist and bList as they are objects
        // ==================================================================
        List<List> AList = new ArrayList<>();
        List<Double> bList = new ArrayList<>();
        boolean status = preProcess(A,b,AList,bList);
        if(!status) return -1;
        
        // ==================== make subsets ================================
        // take the size of Alist and traverse all subsets of sets
        // ==================================================================
        List<List> sets = setTraverse(bList.size(),AList.get(0).size());
        
        System.out.print(Arrays.toString(sets.toArray()));
        int index = 0;
        // ============== rebuild new Equations =============================
        // tak list AList, and bList and one subset of sets
        // generate Equation that contains matrix a and vector b
        // ==================================================================
        List<Integer> subset = sets.get(index);
        Equation myEQ = buildEquations(AList, bList, subset);
        
      
        // ==================== Solve for Vertex "p" ========================
        // ==================================================================
        double[] soln = solveEquations(myEQ);
        System.arraycopy(soln, 0, x, 0, soln.length);
        
        
        
        // ==================== reserved section ================================
        // ==================================================================
        
        // ==================== reserved section ================================
        // ==================================================================
        
      return 0;
    }
    
    static boolean preProcess(double[][] A, double[] b, List<List> AList, List<Double> bList){
        // ==================== Build Initial A,b Lists =====================
        // Matrix (Arrays) A and b can not be modified, transfer these to 
        // arrayLists first
        // ==================================================================
        for(int i =0; i<A.length; i++){
            List RowList = new ArrayList<Double>();
            for(int j=0; j<A[0].length; j++){
                RowList.add(A[i][j]);
            }
            AList.add(new ArrayList(RowList));
            bList.add(b[i]);
        }

//        System.out.println("initial");
//        System.out.println("the matix A: "+ Arrays.toString(AList.toArray()));
//        System.out.println("the vector b: "+ Arrays.toString(bList.toArray()));
//        //System.out.println("the subset chosen: " + Arrays.toString(sets.toArray()));
        
        
        
        // =============== add non-negativity constraints ===================
        // all x's >= 0   this should be written as    -x_i <= 0
        // ==================================================================
        int n_x = AList.get(0).size();
        for(int i=0; i<n_x; i++){
            List newRow = new ArrayList<Double>();
            for(int j=0; j<n_x; j++){
                double element = (j==i) ? -1.0: 0.0;
                newRow.add(element);
            }
            AList.add(new ArrayList(newRow));
            bList.add((double) 0);
        }
        
//        System.out.println("after adding non-negativity");
//        System.out.println("the matix A: "+ Arrays.toString(AList.toArray()));
//        System.out.println("the vector b: "+ Arrays.toString(bList.toArray()));
//        //System.out.println("the subset chosen: " + Arrays.toString(sets.toArray()));
        
        // ==================== Check for Parallel Constraints ==============
        // if A!/A       ==> it is fine, not parallel 
        // if A//A, b//b ==> remove one
        // if A//A, b!/b ==> case 1: one constraint is redundant and should be 
        //                           removed
        //                   case 2: there is no-over lapping of the two 
        //                           ==> report no solution
        //                   case 3: over lapping, that's fine, do nothing
        // ==================================================================
        for(int row =0; row<AList.size(); row++){for(int otherRow = 0; otherRow<AList.size(); otherRow++){if(row != otherRow){
            // ratios is the ratio between each pair of elements of rows in 
            // matrix A 
            //System.out.println("comparing row: "+row+" and other row: "+otherRow);
            List<Double> ratios = new ArrayList<>(); 
            double ARatio = (double) 0.0;
            double bRatio = (double) 0.0;
            
            for(int col =0; col<AList.get(0).size(); col++){
                double tempRatio = (double)AList.get(row).get(col)/(double)AList.get(otherRow).get(col);
                boolean isZero = Objects.equals((double)AList.get(row).get(col), 0.0);
                if (!(Double.isInfinite(tempRatio)||Double.isNaN(tempRatio))||!isZero){
                    ratios.add(tempRatio);
                    //System.out.println(tempRatio);
                }
            }
                bRatio = (double)bList.get(row)/(double)bList.get(otherRow);
            
            if(ratios.stream().distinct().limit(2).count() <= 1){ // A parallel
                //System.out.println("A parallel");
                ARatio = (double) ratios.get(0);
                if(Objects.equals(ARatio, bRatio)){
                    //System.out.println("b parallel");
                    AList.remove(otherRow);
                    bList.remove(otherRow);
                    //System.out.println("*******************removing row: "+otherRow);
                } else if(ARatio > 0){
                    //System.out.println("same direction");
                    if (bList.get(otherRow) > bList.get(row)) {
                        AList.remove(otherRow);
                        bList.remove(otherRow);
                        //System.out.println("*******************removing row: "+otherRow);
                    } else {
                        AList.remove(row);
                        bList.remove(row);
                        //System.out.println("*******************removing row: "+row);
                    }
                } else if ((((double)bList.get(row))+((double)bList.get(otherRow)))<0){
                    //System.out.println("not overlapping");
                    return false; // answer if the system of equations are fine
                }
            }
            
//        System.out.println("after processing");
//        System.out.println("the matix A: "+ Arrays.toString(AList.toArray()));
//        System.out.println("the vector b: "+ Arrays.toString(bList.toArray()));
//        System.out.println();
//        //System.out.println("the subset chosen: " + Arrays.toString(sets.toArray()));
            
            
        }}}
        
        return true;
    }
    
    static List<List> setTraverse(int n, int m){
        List<List> Sets = new ArrayList<>();
        // n # of constraints/restrictions
        // m # of variables
        // m out of n constrains should solve to an point
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
    
    static Equation buildEquations(List<List> AList, List<Double> bList, List<Integer> subset){
        List<ArrayList> AE = new ArrayList(); 
        List<Double> bE = new ArrayList(); 
        
        for(int index =0; index<AList.size(); index++){
            if((int)subset.get(index) == 1){
                AE.add(new ArrayList((AList.get(index))));
                bE.add(bList.get(index));
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
    
    static double[] solveEquations(Equation equation) {
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
    
    void solve() throws IOException {
        int n = nextInt();
        int m = nextInt();
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
        double[] c = new double[m];
        for (int i = 0; i < m; i++) {
            c[i] = nextInt();
        }
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
