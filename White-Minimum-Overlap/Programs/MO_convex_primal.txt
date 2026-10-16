/*********************************************
 * OPL 22.1.0.0 Model
 * Author: ethanwhite
 * Creation Date: Mar 9, 2023 at 6:18:51 PM
 *********************************************/

// The full convex program as recorded in Section 5 of Minimum Overlap Paper https://arxiv.org/abs/2201.05704


 
 int N = 10000;
 int T = 4000;
 int R = 10;
 float L = 2/N;
 
 float h1 = 0;
 float h2 = 0.08;
 float p1 = 0;
 float p2 = 1;
 float q1 = -1;
 float q2 = -0.05;
 
 float pi = 3.1415926535897932384626;
 
 range rangeN = 1..N;
 range rangeT = 1..T;
 range range2T = 1..2*T;
 range rangeR = 1..R;
 range range2R = 1..2*R;
 
 float alphaminus[1..N][1..2*R];
 float alphaplus[1..N][1..2*R];
 float betaplus[1..N][1..2*R];
 float betaminus[1..N][1..2*R];
 
 execute{
   for (var j in rangeN){
     for (var m in range2R){
       alphaminus[j][m] = Math.cos(pi*m*(j-0.5)/N)-pi*m*0.5/N;
       alphaplus[j][m] = Math.cos(pi*m*(j-0.5)/N)+pi*m*0.5/N;
       betaplus[j][m] = Math.sin(pi*m*(j-0.5)/N)+pi*m*0.5/N;
       betaminus[j][m] = Math.sin(pi*m*(j-0.5)/N)-pi*m*0.5/N;
     }
   }
 }
 
 
 
 
 
 float signvec[1..2*T];
 
 execute{ for(var j in range2T){
   signvec[j] = Math.pow(-1,j+1);
 }
   
 }
 
 float epsilonbounds[1..R];
 float deltabounds[1..R];
 
 execute{
   for(var m in rangeR){
     epsilonbounds[m]= (1/pi)*(1/(4-Math.pow(m/T,2)))*2*m*Math.pow(6*T*T*T,-0.5);
     deltabounds[m] = (4/pi)*(1/(4-Math.pow(m/T,2)))*Math.pow(2*T,-0.5);
   }
 }
 
 
 dvar float+ Omega;
 dvar float+ w[1..N];
 dvar float+ v[1..N];
 
 dvar float c[1..T];
 dvar float d[1..T];
 
 dvar float epsilon[1..R];
 dvar float delta[1..R];
 
 dexpr float a_odd[m in 1..R] = epsilon[m]+2*(2*m-1)*signvec[m]*(1/pi)*(0.5/((2*m-1)*(2*m-1))+sum(k in 1..T)(signvec[k+1]*c[k]/((2*m-1)*(2*m-1)-4*k*k)));
 dexpr float a_even[m in 1..R] = 0.5*c[m];
 dexpr float b_odd[m in 1..R] = delta[m]+(4/pi)*signvec[m]*sum(k in 1..T)(k*signvec[k+1]*d[k]/((2*m-1)*(2*m-1)-4*k*k));
 dexpr float b_even[m in 1..R] = 0.5*d[m];
 
 
 minimize Omega;
 
 subject to{
   forall(j in 1..N){
     w[j]<= Omega;
     v[j]<=Omega;
   }
   
   L*sum(j in 1..N)(w[j]+v[j]) == 1;
   
   L*L*sum(j in 1..N)(j*w[j]-(j-1)*v[j])>= h1;
   
   L*L*sum(j in 1..N)(j-1)*(j-1)*(w[j]+v[j]) <= (2/3+h2*h2/2)*N/2;
   
   forall(m in 1..R){
     0.5*L*sum(j in 1..N)alphaminus[j][2*m]*(w[j]+v[j])<= -2*(a_even[m]*a_even[m]+b_even[m]*b_even[m]);
     0.5*L*sum(j in 1..N)alphaminus[j][2*m-1]*(w[j]+v[j])<= 4*signvec[m]*a_odd[m]/(pi*(2*m-1))-2*(a_odd[m]*a_odd[m]+b_odd[m]*b_odd[m]);
     
     0.5*L*sum(j in 1..N)(betaminus[j][2*m]*w[j]-betaplus[j][2*m]*v[j])<= 0;
     0.5*L*sum(j in 1..N)(betaminus[j][2*m-1]*w[j]-betaplus[j][2*m-1]*v[j])<= -signvec[m]*b_odd[m]*8/(pi*(2*m-1));
     
     0.5*L*sum(j in 1..N)(betaplus[j][2*m]*w[j]-betaminus[j][2*m]*v[j])>= 0;
     0.5*L*sum(j in 1..N)(betaplus[j][2*m-1]*w[j]-betaminus[j][2*m-1]*v[j])>= -signvec[m]*b_odd[m]*8/(pi*(2*m-1));
     
     epsilon[m] <= epsilonbounds[m];
     delta[m] <= deltabounds[m];
     epsilon[m] >= -epsilonbounds[m];
     delta[m] >= -deltabounds[m];
   }
   
   forall (k in 1..T){
     c[k]<= 2/pi;
     d[k]<=2/pi;
     c[k]>=-2/pi;
     d[k]>=-2/pi;
   }
   
   sum(k in 1..T)(c[k]*c[k]+d[k]*d[k])<= 1/2;
   
   p1 <= c[1];
   c[1]<= p2;
   q1 <= d[1];
   d[1]<= q2;
   
   0.5*L*sum(j in 1..N)alphaplus[j][2]*(w[j]+v[j])>=-0.5*(p2*p2+q1*q1);
   
 }
 
 
