#include <iostream>
#include <iomanip>
#include <fstream>
#include <dai/alldai.h>
#include <dai/jtree.h>

using namespace std;
using namespace dai;

#define NS_A 3
#define NS_S 5
#define NS_Z 2

#define H 6


// simple variable elimination (eliminate in their default order)
FactorGraph variableElimination(const FactorGraph &fg, const VarSet &elim)
{
  vector<Factor> factors = fg.factors();
  for (VarSet::const_iterator i = elim.begin(); i != elim.end(); i++) {
    Factor f;
    for (vector<Factor>::iterator j = factors.begin(); j!=factors.end(); j++) {
      if (j->vars().contains(*i)) {
	f *= *j;
	j = factors.erase(j)-1;
      }
    }
    factors.push_back(f.marginal(f.vars() / *i));
  }
  return FactorGraph(factors);
}


// add evidence to a Factor Graph
FactorGraph addEvidence(const FactorGraph &fg, const map<Var, size_t> &evid)
{
  vector<Factor> factors = fg.factors();
  for (vector<Factor>::iterator i = factors.begin(); i != factors.end(); i++) {
    for (map<Var,size_t>::const_iterator j = evid.begin(); j!=evid.end(); j++) {
      if (i->vars().contains(j->first)) {
	*i = i->slice(VarSet(j->first), j->second);
      }
    }
  }

  return FactorGraph(factors);
}


// return joint distirbution of a factor graph
Factor jointDistribution(const FactorGraph &fg)
{
  Factor m;
  for (vector<Factor>::const_iterator i = fg.factors().begin();
       i != fg.factors().end(); i++) {
    m *= *i;
  }
  if (m.maxAbs() > 0.0) {
    m.normalize();
  }
  return m;
}


// get var set for all factors in a graph
VarSet allVars(const FactorGraph &fg)
{
  VarSet vs;

  for (vector<Factor>::const_iterator i = fg.factors().begin();
       i != fg.factors().end(); i++) vs |= i->vars();
  return vs;
}

// get maximum value and index from factor
double factorMax(const Factor &f, int *max_i)
{
  double max = f[0];
  *max_i = 0;
  
  for (int i = 1; i < f.states(); i++) {
    if (f[i] > max) {
      max = f[i];
      *max_i = i;
    }
  }
  return max;
}



// value iteration over a factored MDP
void valueIteration(const FactorGraph &dbn, const vector<Factor> &rewards,
		    const VarSet &states, const VarSet &next_states,
		    const VarSet &actions, int horizon, Factor* V, Factor *Pi)
{
  *V = Factor(states, 0.0);
  Factor Q(actions), V2(states);
  
  VarSet elim_vars = ((allVars(dbn) / actions) / states) / next_states;
  FactorGraph r_dbn = variableElimination(dbn, elim_vars);

  for (int t = horizon-1; t >= 0; t--) {
    Pi[t] = Factor(states, 0.0);
    for (State s(states); s.valid(); s++) {
      FactorGraph r_dbn2 = addEvidence(r_dbn, s);
      Q.fill(0.0);
      for (State a(actions); a.valid(); a++) {
	for (int i = 0; i < rewards.size(); i++) {
	  Q[a] += rewards[i][a(rewards[i].vars()) + s(rewards[i].vars())];
	}
	Factor f = jointDistribution(addEvidence(r_dbn2, a));
	Q[a] += (f.p() * V->p()).sum();
      }
      int act;
      V2[s] = factorMax(Q, &act);
      Pi[t][s] = act;
    }
    *V = V2;
  }
}




// compute finite horizon value of a policy
Factor policyValue(const FactorGraph &dbn, const Factor *Pi,
		   const vector<Factor> &rewards, const VarSet &states,
		   const VarSet& next_states, const VarSet &actions,
		   int horizon)
{
  Factor V(states, 0.0), V2(states, 0.0);

  VarSet elim_vars = ((allVars(dbn) / actions) / states) / next_states;
  FactorGraph r_dbn = variableElimination(dbn, elim_vars);

  for (int t = horizon-1; t >= 0; t--) {
    V2.fill(0.0);
    for (State s(states); s.valid(); s++) {
      FactorGraph r_dbn2 = addEvidence(r_dbn, s);
      State a(actions, (int) Pi[t][s]);
      
      for (int i = 0; i < rewards.size(); i++) {
	V2[s] += rewards[i][a(rewards[i].vars()) + s(rewards[i].vars())];
      }
      Factor f = jointDistribution(addEvidence(r_dbn2, a));
      V2[s] += (f.p() * V.p()).sum();
    }
    V = V2;
  }
  return V;
}


// mpc control
Factor mpcValueIteration(const FactorGraph &dbn, const vector<Factor> &rewards,
			 const VarSet *ctrl_states, const VarSet *unctrl_states,
			 const VarSet &actions, const map<Var, size_t> preds,
			 int pred_horizon, int horizon)
{

  // set up DBN with transitions dependent on predictions
  FactorGraph r_dbn = addEvidence(dbn, preds);
  JTree jt(r_dbn, PropertySet("[verbose=0,updates=HUGIN]"));
  jt.init();
  jt.run();

  // set up small DBN (no prediction states)
  VarSet small_vars = actions | ctrl_states[0] | ctrl_states[1] |
    unctrl_states[0] | unctrl_states[1];
  vector<Factor> small_factors;
  for (int i = 0; i < dbn.nrFactors(); i++) {
    if ((dbn.factor(i).vars() & small_vars) == dbn.factor(i).vars()) {
      small_factors.push_back(dbn.factor(i));
    }
  }
  FactorGraph s_dbn(small_factors);

  
  // find all factors that depend entirely on uncontrolled states
  vector<Factor*> s_factors;
  vector<const Factor*> r_factors[pred_horizon];
  for (int i = 0; i < s_dbn.nrFactors(); i++) {
    Factor& f = s_dbn.factor(i);
    if ((f.vars() & (unctrl_states[0] | unctrl_states[1])) == f.vars() &&
	(f.vars() & unctrl_states[0]) != f.vars()) {
      s_factors.push_back(&f);
    }
  }
  for (int t = 0; t < pred_horizon; t++) {
    for (int i = 0; i < dbn.nrFactors(); i++) {
      const Factor &f = dbn.factor(i);
      if ((f.vars() & (unctrl_states[t] | unctrl_states[t+1])) == f.vars() &&
	(f.vars() & unctrl_states[t]) != f.vars()) {
	r_factors[t].push_back(&f);
      }
    }
  }


  // run value iteration
  Factor V(ctrl_states[0] | unctrl_states[0], 0.0);
  Factor V2(ctrl_states[0] | unctrl_states[0], 0.0);
  Factor Pi(ctrl_states[0] | unctrl_states[0], 0.0);
  Factor Q(actions);

  for (int t = horizon-1; t >= 0; t--) {
    
    // use prediction-dependent CPDs for uncontrollable variables
    if (t < pred_horizon) {
      for (int i = 0; i < s_factors.size(); i++) {
	s_factors[i]->p() = jt.calcMarginal(r_factors[t][i]->vars()).p();
      }
    }

    // compute Bellman backup
    for (State s(ctrl_states[0] | unctrl_states[0]); s.valid(); s++) {
      FactorGraph s_dbn2 = addEvidence(s_dbn, s);
      Q.fill(0.0);
      for (State a(actions); a.valid(); a++) {
	for (int i = 0; i < rewards.size(); i++) {
	  Q[a] += rewards[i][a(rewards[i].vars()) + s(rewards[i].vars())];
	}
	Factor f = jointDistribution(addEvidence(s_dbn2, a));
	Q[a] += (f.p() * V.p()).sum();
      }
      int act;
      V2[s] = factorMax(Q, &act);
      Pi[s] = act;
    }
    V = V2;
  }
  return Pi;
}


// get a complete MPC policy for all predictions
Factor mpcPolicy(const FactorGraph &dbn, const vector<Factor> &rewards,
		 const VarSet *ctrl_states, const VarSet *unctrl_states,
		 const VarSet &actions, const VarSet &predictions,
		 int pred_horizon, int horizon)
{
  Factor Pi(ctrl_states[0] | unctrl_states[0] | predictions, 0.0);
  Factor Pi2(ctrl_states[0] | unctrl_states[0]);

  for (State p(predictions); p.valid(); p++) {
    // get MPC policy for this set of predictions
    Pi2 = mpcValueIteration(dbn, rewards, ctrl_states, unctrl_states,
			    actions, p, pred_horizon, horizon);

    for (State s(Pi2.vars()); s.valid(); s++) {
      Pi[p(Pi.vars()) + s(Pi.vars())] = Pi2[s];
    }
  }
  
  return Pi;
}
	

int main(int argc, char* argv[])
{
  Var A, S, Sn, Z[H+2], P[H][H+2];
  int idx = 0;
  
  int horizon = (argc > 1 ? atoi(argv[1]) : 10);
  int pred_horizon = (argc > 2 ? atoi(argv[2]) : H);
  double pred_err = (argc > 3 ? atof(argv[3]) : 0.05);

  // create variables
  A.label() = idx++; A.states() = NS_A;
  S.label() = idx++; S.states() = NS_S;
  Z[0].label() = idx++; Z[0].states() = NS_Z;
  Sn.label() = idx++; Sn.states() = NS_S;
  for (int i = 1; i < H+2; i++) {
    Z[i].label() = idx++; Z[i].states() = NS_Z;
  }
  for (int i = H; i >= 0; i--) {
    for (int j = 0; j < H; j++) {
      if (i + j <= H)
	P[j][i].label() = idx++; P[j][i].states() = NS_Z;
    }
  }

  // create CPDs
  vector<Factor> factors;
  factors.push_back(Factor(VarSet(A, S) | Z[0] | Sn));
  ifstream fin("../code/CPT_S.dat");
  for (int i = 0; i < factors.back().states(); i++) fin >> factors.back()[i];

  
  Prob P_Zn_Z(NS_Z * NS_Z);
  fin.close(); fin.open("../code/CPT_Z.dat");
  for (int i = 0; i < P_Zn_Z.size(); i++) fin >> P_Zn_Z[i];
  
  for (int i = 0; i < H+1; i++) {
    factors.push_back(Factor(VarSet(Z[i], Z[i+1]), P_Zn_Z));
  }

  Prob P_P_Z(NS_Z * NS_Z, pred_err);
  for (int i = 0; i < NS_Z; i++) P_P_Z[i*NS_Z+i] += 1.0 - pred_err*NS_Z;

  for (int i = H; i >= 0; i--) {
    factors.push_back(Factor(VarSet(Z[i+1], P[0][i]), P_P_Z));
    for (int j = 1; j < H; j++) {
      if (j + i <= H) {
	factors.push_back(Factor(VarSet(P[j-1][i+1], P[j][i]), P_P_Z));
      }
    }
  }
  factors.push_back(Factor(A, 1.0/NS_A));
  factors.push_back(Factor(S, 1.0/NS_S));
  factors.push_back(Factor(Z[0], 1.0/NS_Z));
  
  // create DBN, rewards, and state distribution
  FactorGraph dbn(factors);
  vector<Factor> rewards;
  VarSet s(S, Z[0]), sn(Sn, Z[1]);
  Factor V, Pi[horizon], s0, Pi_mpc[horizon], V_mpc;
  s0.normalize();

  for (int i = 0; i < pred_horizon; i++) s |= P[i][0];
  for (int i = 0; i < pred_horizon; i++) sn |= P[i][1];
  rewards.push_back(Factor(VarSet(A,S) | Z[0]));
  fin.close(); fin.open("../code/R.dat");
  for (int i = 0; i < rewards.back().states(); i++) fin >> rewards.back()[i];
  s0 = jointDistribution(variableElimination(dbn, allVars(dbn) / s));

  VarSet ctrl_states[2] = {VarSet(S), VarSet(Sn)};
  VarSet unctrl_states[H+2];
  for (int i = 0; i < pred_horizon+2; i++) unctrl_states[i] = VarSet(Z[i]);
  VarSet predictions;
  for (int i = 0; i < pred_horizon; i++) predictions |= P[i][0];

  valueIteration(dbn, rewards, s, sn, VarSet(A), horizon, &V, Pi);
  cout << setprecision(10) << (V * s0).sum() << endl;

  for (int t = 0; t < horizon; t++) {
    Pi_mpc[t] = mpcPolicy(dbn, rewards, ctrl_states, unctrl_states, VarSet(A),
			  predictions, pred_horizon, horizon-t);
  }
  V_mpc = policyValue(dbn, Pi_mpc, rewards, s, sn, VarSet(A), horizon);
  cout << setprecision(10) << (V_mpc * s0).sum() << endl;

  /*
  cout << V << endl;
  Factor V2(s);
  fin.close(); fin.open("../code/V.dat");
  for (int i = 0; i < V2.states(); i++) fin >> V2[i];
  for (int i = 0; i < V.states(); i++) cout << (V[i]-V2[i])/fabs(V[i]) << endl;
  */

  
  
  

  return 0;
}
  

  
  
  


/*

  
  
  
// simple variable elimination
Factor variableElimination(FactorGraph &fg, const VarSet &query,
			 const map<Var,size_t> &evid)
{
  vector<Var> elim_vars = fg.vars();
  vector<Factor> factors = fg.factors();
  Factor m;

  // find variables to eliminate
  for (vector<Var>::iterator i = elim_vars.begin(); i != elim_vars.end(); i++) {
    if (query.contains(*i)) i = elim_vars.erase(i)-1;
    if (evid.count(*i) > 0) i = elim_vars.erase(i)-1;
  }

  // eliminate variables
  for (vector<Var>::iterator i = elim_vars.begin(); i != elim_vars.end(); i++) {
    Factor f;
    for (vector<Factor>::iterator j = factors.begin(); j!=factors.end(); j++) {
      if (j->vars().contains(*i)) {
	f *= *j;
	j = factors.erase(j)-1;
      }
    }
    factors.push_back(f.marginal(f.vars() / *i));
  }

  // add evidence
  for (vector<Factor>::iterator i = factors.begin(); i != factors.end(); i++) {
    for (map<Var,size_t>::const_iterator j = evid.begin(); j!=evid.end(); j++) {
      if (i->vars().contains(j->first)) {
	*i = i->slice(VarSet(j->first), j->second);
      }
    }
  }

  // multiply remaining marginals
  for (vector<Factor>::iterator i = factors.begin(); i != factors.end(); i++) {
    m *= *i;
  }
  m.normalize();
  return m;
}
*/


 /*
  JTree jt(dbn, PropertySet("[verbose=1,updates=HUGIN]"));
  VarSet s(A, S); s |= Z[0]; for (int i = 0; i < H; i++) s |= P[i][0];
  VarSet sn(Sn,Z[1]); for (int i = 0; i < H; i++) sn |= P[i][1];
  
  map<Var,size_t> evid;
  evid[A] = 0;
  evid[S] = 0;
  evid[Z[0]] = 0;
  for (int i = 0; i < H; i++) evid[P[i][0]] = 0;
  

  for (map<Var,size_t>::const_iterator i = evid.begin(); i != evid.end(); i++)
    jt.clamp(dbn.findVar(i->first), i->second);
  jt.init();
p  jt.run();
  
  FactorGraph rdbn = variableElimination(dbn, (allVars(dbn)/s)/sn);
  Factor f3 = jointDistribution(addEvidence(rdbn, evid));
 */
