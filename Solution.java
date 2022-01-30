public class Solution {
	
	static long gcd(long a, long b)
	{
	    if (a == 0)
	        return b;
	    else if (b == 0)
	        return a;
	    if (a < b)
	        return gcd(a, b % a);
	    else
	        return gcd(b, a % b);
	}
	     
	// Function to convert decimal to fraction
	static String decimalToFraction(double number)
	{
	   
	    // Fetch integral value of the decimal
	    double intVal = Math.floor(number);
	    System.out.println("NUMBER:"+number);
	    // Fetch fractional part of the decimal
	    double fVal = number - intVal;
	    int checker=100;
	    double frac;
	    if(number>1)
	    {
	    	int real=(int)Math.round(number);
	    	frac=number-real;
	    }
	    else
	    	frac=number;
	    	for(int i=2;i<200;i++)
	    	{
	    		if(Math.round(i*frac)-i*frac<=0.005)
	    		if(Math.abs((i*frac)-Math.round(i*frac))<0.02)
	    		{
	    			checker=i;
	    	//		System.out.print("CHECKER:"+checker);
	    			break;
	    		}
	    	}
	    	if(checker>100)
	    		checker=100;
	    	final long pVal=checker;
	    // Consider precision value to
	    // convert fractional part to
	    // integral equivalent  
	    // Calculate GCD of integral
	    // equivalent of fractional
	    // part and precision value
	    long gcdVal = gcd(Math.round(
	                      fVal * pVal), pVal);
	   
	    // Calculate num and deno
	    long num = Math.round(fVal * pVal) / gcdVal;
	    long deno = pVal / gcdVal;
	    return (long)(intVal * deno) +
                num + "/" + deno;
	}
	  public static int LCM(int[] element_array)
	    {
		  int[] element_array1=element_array;
	        int lcm_of_array_elements = 1;
	        int divisor = 2;
	         
	        while (true) {
	            int counter = 0;
	            boolean divisible = false;
	             
	            for (int i = 0; i < element_array1.length; i++) {
	 
	                // lcm_of_array_elements (n1, n2, ... 0) = 0.
	                // For negative number we convert into
	                // positive and calculate lcm_of_array_elements.
	 
	                if (element_array1[i] == 0) {
	                    return 0;
	                }
	                else if (element_array1[i] < 0) {
	                    element_array1[i] = element_array1[i] * (-1);
	                }
	                if (element_array1[i] == 1) {
	                    counter++;
	                }
	 
	                // Divide element_array by devisor if complete
	                // division i.e. without remainder then replace
	                // number with quotient; used for find next factor
	                if (element_array1[i] % divisor == 0) {
	                    divisible = true;
	                    element_array1[i] = element_array1[i] / divisor;
	                }
	            }
	 
	            // If divisor able to completely divide any number
	            // from array multiply with lcm_of_array_elements
	            // and store into lcm_of_array_elements and continue
	            // to same divisor for next factor finding.
	            // else increment divisor
	            if (divisible) {
	                lcm_of_array_elements = lcm_of_array_elements * divisor;
	            }
	            else {
	                divisor++;
	            }
	 
	            // Check if all element_array is 1 indicate
	            // we found all factors and terminate while loop.
	            if (counter == element_array1.length) {
	            	if(lcm_of_array_elements<100)
	                return lcm_of_array_elements;
	            	else
	            	{
	            		return 100;
	            	}
	            }
	        }
	    }
	
	public static double[][] MultiplyMatrices(double[][] A, double[][] B) {

        int aRows = A.length;
        int aColumns = A[0].length;
        int bRows = B.length;
        int bColumns = B[0].length;
        double[][] C = new double[aRows][bColumns];
        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                C[i][j] = 0.00000;
            }
        }

        for (int i = 0; i < aRows; i++) { // aRow
            for (int j = 0; j < bColumns; j++) { // bColumn
                for (int k = 0; k < aColumns; k++) { // aColumn
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }
	
	
	 public static double[][] invert(double a[][]) 
	    {
	        int n = a.length;
	        double x[][] = new double[n][n];
	        double b[][] = new double[n][n];
	        int index[] = new int[n];
	        for (int i=0; i<n; ++i) 
	            b[i][i] = 1;
	 
	 // Transform the matrix into an upper triangle
	        gaussian(a, index);
	 
	 // Update the matrix b[i][j] with the ratios stored
	        for (int i=0; i<n-1; ++i)
	            for (int j=i+1; j<n; ++j)
	                for (int k=0; k<n; ++k)
	                    b[index[j]][k]
	                    	    -= a[index[j]][i]*b[index[i]][k];
	 
	 // Perform backward substitutions
	        for (int i=0; i<n; ++i) 
	        {
	            x[n-1][i] = b[index[n-1]][i]/a[index[n-1]][n-1];
	            for (int j=n-2; j>=0; --j) 
	            {
	                x[j][i] = b[index[j]][i];
	                for (int k=j+1; k<n; ++k) 
	                {
	                    x[j][i] -= a[index[j]][k]*x[k][i];
	                }
	                x[j][i] /= a[index[j]][j];
	            }
	        }
	        return x;
	    }
	    public static void gaussian(double a[][], int index[]) 
	    {
	        int n = index.length;
	        double c[] = new double[n];
	 
	 // Initialize the index
	        for (int i=0; i<n; ++i) 
	            index[i] = i;
	 
	 // Find the rescaling factors, one from each row
	        for (int i=0; i<n; ++i) 
	        {
	            double c1 = 0;
	            for (int j=0; j<n; ++j) 
	            {
	                double c0 = Math.abs(a[i][j]);
	                if (c0 > c1) c1 = c0;
	            }
	            c[i] = c1;
	        }
	        int k = 0;
	        for (int j=0; j<n-1; ++j) 
	        {
	            double pi1 = 0;
	            for (int i=j; i<n; ++i) 
	            {
	                double pi0 = Math.abs(a[index[i]][j]);
	                pi0 /= c[index[i]];
	                if (pi0 > pi1) 
	                {
	                    pi1 = pi0;
	                    k = i;
	                }
	            }
	            int itmp = index[j];
	            index[j] = index[k];
	            index[k] = itmp;
	            for (int i=j+1; i<n; ++i) 	
	            {
	                double pj = a[index[i]][j]/a[index[j]][j];
	                a[index[i]][j] = pj;
	                for (int l=j+1; l<n; ++l)
	                    a[index[i]][l] -= pj*a[index[j]][l];
	            }
	        }
	        }

	 public static int[] solution(int[][] m) {
		 if(m.length==1)
		 {
			 int arr[]= {1,1};
			 return arr;
		 }
		 int[][] nullMatrix = new int[m.length][m[0].length];
	        int[][] notNullMatrix=new int[m.length][m[0].length];	//USING MARKOV CHAINS 
	        int arrayOfNullIndexes[]=new int[m.length];
	        int ctrNull=0;
	        int ctrNotNull=0;
		 for(int i=0;i<m.length;i++)
		 {
			 boolean checkIfNull=true;
			 for(int j=0;j<m[0].length;j++)
			 {
				 if(m[i][j]!=0)
					 checkIfNull=false;
			 }
			 if(checkIfNull==true)
			 {
				 nullMatrix[ctrNull++]=m[i];
				 arrayOfNullIndexes[i]=1;
			 }
			 else
			 {
				 notNullMatrix[ctrNotNull++]=m[i];
				 arrayOfNullIndexes[i]=0;
			 }
		 }
		 int[][] finalNullMatrix=new int[ctrNull][ctrNull];
		 int[][] finalNotNullMatrix=new int[ctrNotNull][ctrNotNull];
	        int[] denominators=new int[ctrNotNull];
	        
	        
		 for(int i=0;i<ctrNull;i++)			//FINAL CLEAN MATRICES FOR NULL AND NOT NULL, NULL MATRIX BEING THE IDENTITY MATRIX
		 {
			 finalNullMatrix[i]=nullMatrix[i];			//IZDVOJI MATRICU GDE SE NALAZE STANJA SAMO ZA NOTNULL
		 }
		 for(int i=0;i<ctrNotNull;i++)
		 {
			 finalNotNullMatrix[i]=notNullMatrix[i];
		 }
		 for(int i=0;i<ctrNull;i++)
		 {
			 finalNullMatrix[i][i]=1;
		 }
		 for(int i=0;i<ctrNotNull;i++)
		 {
			 int sum=0;
			 for(int j=0;j<m.length;j++)
			 {
				 sum+=finalNotNullMatrix[i][j];
			 }
			 denominators[i]=sum;		//setting the Denominators 
		 }
					//END CHECK
		 int[][] QMatrix=new int[ctrNotNull][ctrNotNull];
		 int qXCounter=0;
		 int qYCounter=0;
		 for(int i=0;i<ctrNotNull;i++)
		 {
			 for(int j=0;j<m.length;j++)
			 {
				if(arrayOfNullIndexes[j]==0)
				{
					QMatrix[qXCounter][qYCounter++]=finalNotNullMatrix[i][j];
				}
			 }
			 qXCounter++;
			 qYCounter=0;
		 }
		 int[][] identityMatrix=new int[ctrNotNull][ctrNotNull];
		 for(int i=0;i<ctrNotNull;i++)
		 {
			 for(int j=0;j<ctrNotNull;j++)
			 {
				 identityMatrix[i][i]=1;
			 }
		 }
		 		//NO NEED FOR NULL MATRIX FROM THE START
		 //CALCULATE I-Q
		 double[][] IQM=new double[ctrNotNull][ctrNotNull];
		 for(int i=0;i<ctrNotNull;i++)
		 {
			 for(int j=0;j<ctrNotNull;j++)
			 {
				 double d= new Double(QMatrix[i][j]);
				 if(QMatrix[i][j]!=0)
				 IQM[i][j]=identityMatrix[i][j]-d/denominators[i];
				 else
					 IQM[i][j]=identityMatrix[i][j];		//IQM CAST TO DOUBLE TO SAVE THE FRACTIONS
			 }
		 }
		IQM= Solution.invert(IQM);
		 double[][] RMatrix=new double[ctrNotNull][ctrNull];
		 int rXCounter=0;
		 int rYCounter=0;
		 for(int i=0;i<ctrNotNull;i++)
		 {
			 for(int j=0;j<m.length;j++)
			 {
				if(arrayOfNullIndexes[j]==1)
				{
					double d=new Double(denominators[i]);
					RMatrix[rXCounter][rYCounter++]=finalNotNullMatrix[i][j]/d;		//R MATRIX WORKING
				}
			 }
			 rXCounter++;
			 rYCounter=0;
		 }
		 
		 
		 double Final[][]=Solution.MultiplyMatrices(IQM, RMatrix);
		 String[] stringArr=new String[Final[0].length];
		 int[] denoms=new int[Final[0].length];
		 int[] noms=new int[Final[0].length];
		 for(int i=0;i<Final[0].length;i++)
		 {
			 stringArr[i]=decimalToFraction(Final[0][i]);
			 System.out.print(" "+stringArr[i]);
		 }
		for(int i=0;i<stringArr.length;i++)
		{
			String[] temp=stringArr[i].split("/");
			denoms[i]=Integer.parseInt(temp[1]);
			noms[i]=Integer.parseInt(temp[0]);
		}
		int LCMnum=Solution.LCM(denoms);
		System.out.println(LCMnum);
			for(int i=0;i<stringArr.length;i++)
			{
				String[] temp=stringArr[i].split("/");	//HAD TO MAKE DENUMS AGAIN
				denoms[i]=Integer.parseInt(temp[1]);
			}
		int[] lastArr=new int[noms.length+1];
		for(int i=0;i<lastArr.length-1;i++)
		{
			double d= (double) LCMnum/denoms[i];
			lastArr[i]=(int) Math.round(noms[i]*d);
		}
		lastArr[noms.length]=LCMnum;
		 for(int i=0;i<lastArr.length;i++)
			 System.out.print(" "+lastArr[i]);
		 return lastArr;
	    }
	public static void main(String[] args) {
		int[][] f= {
				 {0, 86, 61, 189, 0, 18, 12, 33, 66, 39},
			        {0, 0, 2, 0, 0, 1, 0, 0, 0, 0},
			        {15, 187, 0, 0, 18, 23, 0, 0, 0, 0},
			        {1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
			        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		};
		int n[]=Solution.solution(f);

	

}
}
