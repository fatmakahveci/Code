public class FindElementInSortedMatrix {
    public static void main(String[] args) {
        int[][] m = new int[4][4];
        m[0][0]=10;m[0][1]=20;m[0][2]=30;m[0][3]=40;
        m[1][0]=15;m[1][1]=25;m[1][2]=35;m[1][3]=45;
        m[2][0]=27;m[2][1]=29;m[2][2]=37;m[2][3]=48;
        m[3][0]=32;m[3][1]=33;m[3][2]=39;m[3][3]=50;

        int element=25;

        int i=0, j=3;
        while(i>=0 && j<=4) {
            if(element==m[i][j])
                break;
            else if(element<m[i][j])
                j--;
            else
                i++;
        }
        System.out.println("i: "+i+" "+"j: "+j);
    }
}
