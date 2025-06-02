package paper1.NE0602.NE4;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;

public class mnCalForCCSBIZ {
    private int k;
    private int d;
    private double alpha;
    private ArrayList<Double> gammas;
    private ArrayList<Double> izsp = new ArrayList<Double>();  // feasibility tolerance parameter
    private int n0;
    private ArrayList<Integer> mbetan0 = new ArrayList<Integer>();
    private double alpha2;

    public mnCalForCCSBIZ(double alpha, int d, int k, ArrayList<Double> gammas, ArrayList<Double> izsp){
        this.alpha = alpha;
        this.k = k;
        this.d = d;
        this.gammas = gammas;
        this.izsp.addAll(izsp);
    }

    public int getn0(){return n0;}
    public ArrayList<Integer> getmbetan0(){return mbetan0;}
    public double getAlpha2(){return alpha2;}

    public void run(){
        double[] as =  new double[d];
        for(int i = 0 ; i < d ; i++) {
            as[i] = (Math.sqrt((gammas.get(i)-izsp.get(i))*(1-gammas.get(i)+izsp.get(i)))+Math.sqrt(gammas.get(i)*(1-gammas.get(i))))/(izsp.get(i));
        }
//        System.out.println(Arrays.toString(as));
        double[] beta1s = beta1sCalculate(as, alpha, d);
//        System.out.println(Arrays.toString(beta1s));
        alpha2 = beta1s[0];
        for (int i = 1; i < d; i++){
            if (beta1s[i]>alpha2){
                alpha2 = beta1s[i];
            }
        }
//        System.out.println("        "+(alpha2*(k-1)+beta1s[0]+beta1s[1]+beta1s[2]));
//        System.out.println("        "+alpha2*(k-1)+" "+beta1s[0]+" "+beta1s[1]+" "+beta1s[2]);
//        System.out.println(alpha2+" "+alpha/k);

        double[] nsbeta1s =  new double[d];
        for (int i = 0; i < d; i++) {
            nsbeta1s[i] = inverseCumulativeProbabilityforSTD(1 - beta1s[i]) * inverseCumulativeProbabilityforSTD(1 - beta1s[i]) / (izsp.get(i) * izsp.get(i))
                    * (Math.sqrt((gammas.get(i) - izsp.get(i)) * (1 - gammas.get(i) + izsp.get(i))) + Math.sqrt(gammas.get(i) * (1 - gammas.get(i))))
                    * (Math.sqrt((gammas.get(i) - izsp.get(i)) * (1 - gammas.get(i) + izsp.get(i))) + Math.sqrt(gammas.get(i) * (1 - gammas.get(i))));
        }

        double nsbeta1smax = nsbeta1s[0];
        for (int i = 1; i < d; i++){
            if (nsbeta1s[i]>nsbeta1smax){
                nsbeta1smax = nsbeta1s[i];
            }
        }
        n0 = (int)Math.ceil(nsbeta1smax);
//        int[] mbetan0 =  new int[d];
        double mbetan01;
        for(int i = 0 ; i < d ; i++) {
            mbetan01 = n0*gammas.get(i) - inverseCumulativeProbabilityforSTD(1-beta1s[i])*Math.sqrt(n0*gammas.get(i)*(1-gammas.get(i)));
            mbetan0.add((int)Math.floor(mbetan01));
        }
    }

    public double[] beta1sCalculate(double[] as, double alpha, int d){
        double asbar = as[0];
        double asmin = as[0];
        for (int i = 1; i < d; i++){
            if (as[i]>asbar){
                asbar = as[i];
            }
            if (as[i]<asmin){
                asmin = as[i];
            }
        }
        double xi = 10000000;
//        double lastvalue= 10000000;

        double left = -asmin*cumulativeProbabilityforSTD(1); //+inf
        System.out.println(cumulativeProbabilityforSTD(1));
        double right = -asbar/2;  //-inf
//        System.out.println("left: "+left+" right: "+right);
        double[] beta1s =  new double[d];
        double aaright = (k-1)*inverseCumulativeProbabilityforSTD(-right/asbar) - alpha;
        double aaleft = (k-1)*inverseCumulativeProbabilityforSTD(-left/asbar) - alpha;
        for (int j = 0; j < d; j++){
            aaright = aaright + inverseCumulativeProbabilityforSTD(-right/as[j]);
            aaleft = aaleft + inverseCumulativeProbabilityforSTD(-left/as[j]);
        }
        if (left > right || aaleft < 0 || aaright > 0){
            System.out.println("cant calculate n0");
            return beta1s;
        }
        double solutioninterval = right - left;
        int bscount = 0;
//        System.out.println(bscount+"left: "+left+" right: "+right+" solutioninterval: "+solutioninterval);
        while (solutioninterval > 0.000000000001){

            bscount =+1;
            double mid = (left + right)/2;
            double a = (k-1)*inverseCumulativeProbabilityforSTD(-mid/asbar) - alpha;
            for (int j = 0; j < d; j++){
                a = a + inverseCumulativeProbabilityforSTD(-mid/as[j]);
            }
            if (a==0){
//                xi = mid;
                break;
            }else if (a<0){
                right = mid;
            }else{
                left = mid;
            }
            solutioninterval = right - left;
//            System.out.println(bscount+"left: "+left+" right: "+right+" solutioninterval: "+solutioninterval);
        }
        xi = (left + right)/2;
//        double[] beta1s =  new double[d];
//        double[] nsbeta1s =  new double[d];
        for (int i = 0; i < d; i++){
            beta1s[i] = inverseCumulativeProbabilityforSTD(-xi/as[i]);
//            System.out.println(-xi/as[i]);
//            System.out.println(inverseCumulativeProbabilityforSTD(-xi/as[i]));
//            System.out.println(beta1s[i]);
//            nsbeta1s[i] = inverseCumulativeProbabilityforSTD(1-beta1s[i])*inverseCumulativeProbabilityforSTD(1-beta1s[i])/(izsp.get(i)*izsp.get(i))
//                    *(Math.sqrt((gammas.get(i)-izsp.get(i))*(1-gammas.get(i)+izsp.get(i)))+Math.sqrt(gammas.get(i)*(1-gammas.get(i))))
//                    *(Math.sqrt((gammas.get(i)-izsp.get(i))*(1-gammas.get(i)+izsp.get(i)))+Math.sqrt(gammas.get(i)*(1-gammas.get(i))));
        }

        return beta1s;
//        double nsbeta1smax = nsbeta1s[0];
//        for (int i = 1; i < d; i++){
//            if (nsbeta1s[i]>nsbeta1smax){
//                nsbeta1smax = nsbeta1s[i];
//            }
//        }
//        return  nsbeta1smax;
    }

    public double inverseCumulativeProbabilityforSTD(double p){
        NormalDistribution normalDistribution = new NormalDistribution();
        return normalDistribution.inverseCumulativeProbability(p);
    }

    public double cumulativeProbabilityforSTD(double x){
        NormalDistribution normalDistribution = new NormalDistribution();
        return normalDistribution.cumulativeProbability(x);
    }
}
