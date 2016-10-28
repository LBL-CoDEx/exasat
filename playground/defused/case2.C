void func(float *E,float *Eprev,float *R,int n,double alpha,double dt)
{
  int _idx;
  int _gidx;
  _idx = n;
  _gidx = _idx + n;
  int _idy;
  int _gidy;
  _idy = n;
  _gidy = _idy + n;
  {
    if (_gidy >= 1 && _gidy <= n) {{
        if (_gidx >= 1 && _gidx <= n) {
 
	  float m; 

	  E[_gidx] = Eprev[_gidx] + alpha;

	  m=1 ;
          (E[_gidy]) = m + (-dt * (((((((E[_gidy])[_gidx)) * (((E[_gidy])[_gidx]) - dt)) * (((E[_gidy])[_gidx]) - (1))) + (((E[_gidy])[_gidx]) * ((R[_gidy])[_gidx])))));
	  // R[j][i] += dt*(epsilon+M1* R[j][i]/( E[j][i]+M2))*(-R[j][i]-kk*E[j][i]*(E[j][i]-b-1));


        }
      }
    }
  }
}
