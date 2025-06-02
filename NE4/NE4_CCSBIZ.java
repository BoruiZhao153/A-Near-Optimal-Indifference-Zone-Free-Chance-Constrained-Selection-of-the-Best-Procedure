package paper1.NE0602.NE4;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.ArrayList;
import java.util.Random;

public class NE4_CCSBIZ {
    private int k;
    private int d;
    private double alpha;
    private ArrayList<Double> gammas;
    private double izpp;  //indifference-zone parameter
    private ArrayList<Double> izsp = new ArrayList<Double>();  // feasibility tolerance parameter
    private double feasibleThreshold;
    private int bestID = -1;
    private double totalSampleSize=0;
    private ArrayList<int[]> listALT;
    private Random R= new Random();
    private int n0;
    private ArrayList<Integer> mbetan0 = new ArrayList<Integer>();
    private double alpha2;


    public NE4_CCSBIZ(double alpha, int d, ArrayList<Double> gammas, double feasibleThreshold
            , ArrayList<int[]> listALT, long seed, double izpp
            , ArrayList<Double> izsp, int n0, ArrayList<Integer> mbetan0, double alpha2){
        this.alpha = alpha;
        this.d = d;
        this.gammas = gammas;
        this.feasibleThreshold = feasibleThreshold;
        this.listALT = listALT;
        R.setSeed(seed);
        this.izpp = izpp;
        this.izsp.addAll(izsp);
        this.n0 = n0;
        this.mbetan0 = mbetan0;
        this.alpha2 = alpha2;
    }
    public NE4_CCSBIZ(){
    }

    public int getBestID() {
        return bestID;
    }

    public double getTotalSampleSize() {
        return totalSampleSize;
    }

    public void run() {
        k = listALT.size();
        NE4_ED simulation = new NE4_ED(R.nextLong());
//        System.out.println(alpha2);
        long startTime1 = System.currentTimeMillis();
        ArrayList<Integer> I = new ArrayList<Integer>();
        for(int i = 0 ; i < k ; i++) {
            I.add(i);
        }

        double h2 = (n0 - 1)*(Math.pow(2*alpha2,-2.0/(n0-1))-1);

        // Feasibility check
        int fcellnum = 0;
        int ocellnum = 0;
        ArrayList<Integer> F = new ArrayList<Integer>(); //I feasible
        ArrayList<Integer> Fold = new ArrayList<Integer>(); //I feasible
        String[] altFeasible = new String[k];
        for (int i = 0; i < altFeasible.length; i++) {
            altFeasible[i] = "unknown";
        }
        int i1 = 0;
        double[][] sampleinfea =  new double[k][n0];
        while (i1<k){
            int tau = 0;
            int[] Ztau = new int[d];
            while (altFeasible[i1].equals("unknown")){
                tau = tau + 1;
                simulation.setAlt(listALT.get(I.get(i1)));
                simulation.run();
                sampleinfea[i1][tau-1]= -simulation.getAverageWTForC4C5();
//                double samplesp = R.nextGaussian() * Math.sqrt(spsigma2.get(I.get(i1)))+spmu.get(I.get(i1))-feasibleThreshold; /////////////////////////////
                boolean anyViolate = false;
//                    if (i1==0){
//                        System.out.println(sampleinfea[i1][tau-1]);
//                    }
                for (int j = 0; j < d; j++) {
                    double sampleSPp= 0.0;
                    if (j==0){
                        sampleSPp= simulation.getSatisfyC1OrNot();
                    }else if (j==1){
                        sampleSPp= simulation.getSatisfyC2OrNot();
                    }else {
                        sampleSPp= simulation.getSatisfyC3OrNot();
                    }

                    if(sampleSPp<feasibleThreshold){
                        Ztau[j] = Ztau[j] + 1;
                    }
                    if (Ztau[j] >= mbetan0.get(j)+1){
                        anyViolate = true;
//                        System.out.println(tau+"时方案"+I.get(i1)+"因约束"+j+"被判定为不可行"+" 其Ztau均值为："+(double)Ztau[j]/tau+" Ztau: "+Ztau[j]+" tau: "+ mbetan0.get(j));
                    }
                }


                if(anyViolate){
//                    System.out.println("Ztau: "+Ztau+" tau: "+tau);
                    altFeasible[i1] = "infeasible";
//                    fcellnum = fcellnum + 1;
//                    System.out.println(tau+"时长通过可行性分析淘汰了"+I.get(i1)+" 可行性分析已淘汰" + fcellnum);
                    totalSampleSize = totalSampleSize + tau;
//                    numsamplefeasibleell = numsamplefeasibleell + tau;
                }else if(tau == n0){
                    ////////////////////////////////////////////////////////////////////////
                    altFeasible[i1] = "feasible";
//                    numsamplefeasiblepass = numsamplefeasiblepass+tau;
//                    System.out.println(tau+"时长"+I.get(i1)+"通过了可行性分析");
                    F.add(i1);
                }
            }
            i1 = i1 + 1;
        }
//        System.out.println("可行方案集F的大小为： "+F.size()+" 可行性分析进行淘汰花费的样本数为："+numsamplefeasibleell+" 通过可行性分析花费的样本数为："+numsamplefeasiblepass);
        long endTime2 = System.currentTimeMillis();
        long runTime2 = endTime2 - startTime1;
//        System.out.println("可行性分析后方案集F的大小为： "+F.size()+" 样本量为："+totalSampleSize+" 花费时间为："+runTime2);

        //Selecting the Best
        double[] tempX = new double[F.size()];
        double[][] Xij =  new double[k][k];
        double[][] Xij2 = new double[k][k];
        double[][] Sij2 = new double[k][k];
        for(int ell = 0; ell < n0; ell++) {
            for(int i = 0; i < F.size(); i++) {
                tempX[i]= sampleinfea[F.get(i)][ell];
            }

            for(int i = 0; i < F.size(); i++) {
                Xij[F.get(i)][F.get(i)] = Xij[F.get(i)][F.get(i)] + tempX[i];
                Xij2[F.get(i)][F.get(i)] = Xij2[F.get(i)][F.get(i)] + tempX[i] * tempX[i];

                for(int j = i + 1; j < F.size(); j++) {
                    Xij[F.get(i)][F.get(j)] = Xij[F.get(i)][F.get(j)]+tempX[i]-tempX[j];
                    Xij[F.get(j)][F.get(i)] = Xij[F.get(j)][F.get(i)]+tempX[j]-tempX[i];
                    Xij2[F.get(i)][F.get(j)] = Xij2[F.get(i)][F.get(j)] +(tempX[i]-tempX[j])*(tempX[i]-tempX[j]);
                    Xij2[F.get(j)][F.get(i)] = Xij2[F.get(i)][F.get(j)];
                }
            }
        }

        for(int i = 0; i < F.size(); i++) {
            for(int j = i + 1; j < F.size(); j++) {
                Sij2[F.get(i)][F.get(j)] = (Xij2[F.get(i)][F.get(j)] - Xij[F.get(i)][F.get(j)]* Xij[F.get(i)][F.get(j)]/n0)/(n0-1);
                Sij2[F.get(j)][F.get(i)] = Sij2[F.get(i)][F.get(j)];
            }
        }

        int r = n0;

        while(F.size()>1) {
            for(int i  = 0;  i < F.size(); i++) {
                for (int j = 0; j < F.size(); j++) {
                    if(i!=j) {
                        double Zij = Xij[F.get(i)][F.get(j)];
                        double gt = h2*Sij2[F.get(i)][F.get(j)]/(2.0*izpp)-izpp*r/2.0;
                        if (gt < 0){
//                            System.out.println("gt<0000000000000000000000000000000000000000000000000000000000000000000000");
                            gt = 0;
                        }
                        if(Zij > gt) {
                            totalSampleSize = totalSampleSize + r;
//                            numsampleoptimalonly = numsampleoptimalonly + r - n0;
//                            ocellnum = ocellnum + 1;
//                            if (r == n0){
//                                numaltinoptimalwaste = numaltinoptimalwaste + 1;
//                            }
//                            System.out.println(r+"时"+F.get(i)+"凭借pp为"+Xij[F.get(i)][F.get(i)]/r+"的表现"+"通过最优性分析去掉了"+F.get(j)+" 已经通过可行性分析淘汰"+fcellnum
//                                    +" 已通过最优性分析淘汰"+ocellnum+" 总共淘汰了"+(fcellnum+ocellnum)+" 其pp为："+Xij[F.get(j)][F.get(j)]/r);

                            F.remove(j);
                            if(j < i) {
                                i--;
                            }
                            j--;
                        }
                    }
                }
            }
            if(F.size()==1) {
                totalSampleSize = totalSampleSize + r;
//                numsampleoptimalonly = numsampleoptimalonly + r - n0;
            }
            r = r+1;
            double[] tempX2 = new double[F.size()];
            for(int i = 0; i < F.size(); i++) {
                simulation.setAlt(listALT.get(F.get(i)));
                simulation.run();
                tempX2[i]= -simulation.getAverageWTForC4C5();
            }

            for(int i = 0; i < F.size(); i++) {
                Xij[F.get(i)][F.get(i)] = Xij[F.get(i)][F.get(i)] + tempX2[i];
            }

            for(int i = 0; i < F.size(); i++) {
                for(int j = i + 1; j < F.size(); j++) {
                    Xij[F.get(i)][F.get(j)] = Xij[F.get(i)][F.get(j)]+tempX2[i]-tempX2[j];
                    Xij[F.get(j)][F.get(i)] = Xij[F.get(j)][F.get(i)]+tempX2[j]-tempX2[i];
                }
            }
        }
//        System.out.println("可行性分析进行淘汰花费的样本数为："+numsamplefeasibleell+" 通过可行性分析花费的样本数为："+numsamplefeasiblepass
//                +" 最优性分析使用的额外样本数为："+numsampleoptimalonly+" n0对于最优性分析已经过剩的方案数为："+numaltinoptimalwaste+"占此阶段总淘汰方案数的比例为："+(double)numaltinoptimalwaste/ocellnum);
        if(F.size()==1){
            bestID = F.get(0);
        }
        long endTime1 = System.currentTimeMillis();
        long runTime1 = endTime1 - startTime1;
        System.out.println("最优方案为："+bestID+" 样本量为："+totalSampleSize+" 花费时间为："+runTime1);


    }

    public double inverseCumulativeProbabilityforSTD(double p){
        NormalDistribution normalDistribution = new NormalDistribution();
        return normalDistribution.inverseCumulativeProbability(p);
    }

    public double bernoullisample(double p){
//        System.out.println(R.nextDouble());
        boolean result = R.nextDouble() < p;

        if (result) {
            return 1;
        } else {
            return 0;
        }
    }

    private static int factorial(int n) {
        if (n == 0) {
            return 1;
        }
        int result = 1;
        for (int i = 1; i <= n; i++) {
            result *= i;
        }
        return result;
    }

    // 计算组合数 C(n, k) = n! / (k! * (n-k)!)
    private static double binomialCoefficient(int n, int k) {
        if (k > n) return 0;
//        System.out.println(n+" "+k+" "+(n-k));
//        System.out.println((factorial(k)+" "+factorial(n - k)));
//        System.out.println((factorial(k) * factorial(n - k)));
        return (double) factorial(n) / ( factorial(k) * factorial(n - k));
    }

    // 计算二项分布的累积概率 F(x; n, p)
    public static double binomialCDF(int x, int n, double p) {
        double sum = 0.0;
        for (int i = 0; i <= x; i++) {
            sum += binomialCoefficient(n, i) * Math.pow(p, i) * Math.pow(1 - p, n - i);
        }
        return sum;
    }
}
