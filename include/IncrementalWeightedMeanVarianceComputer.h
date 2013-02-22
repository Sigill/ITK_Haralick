#ifndef INCREMENTALWEIGHTEDMEANVARIANCECOMPUTER_H
#define INCREMENTALWEIGHTEDMEANVARIANCECOMPUTER_H

/*
  For weighted distribution, uses West1979, "Updating mean and variance 
  estimates: An improved method".
  SW(1) = w(1), SW(k) = SW(k-1) + SW(k)
  M(1) = x(1), M(k) = M(k-1) + (x(k) - M(k-1) ) * w(k) / SW(k)
  S(1) = 0, S(k) = S(k-1) + (x(k) - M(k-1)) * (x(k) - M(k-1)) * SW(k-1) * W(k) / SW(k)
  for 2 <= k <= n, then
  sigma = vcl_sqrt(S(n) / SW(n))

  This is the version from wikipedia
  http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Weighted_incremental_algorithm
*/
template < typename TValueType >
class IncrementalWeightedMeanVarianceComputer {
public:
	typedef TValueType ValueType;

	IncrementalWeightedMeanVarianceComputer() : SW(0.0), M(0.0), S(0.0) {}
	~IncrementalWeightedMeanVarianceComputer() {}

	void insert(const ValueType value, const ValueType weight)
	{
		if(weight > 0) {
			SW_tmp = SW + weight;
			delta = value - M;
			r = delta * weight / SW_tmp;
			M += r;
			S += SW * delta * r;
			SW = SW_tmp;
		}
	}

	void reset()
	{
		SW = 0.0;
		M = 0.0;
		S = 0.0;
	}

	ValueType getMean() { return M; }
	ValueType getVariance() { return S / SW; }

private:
	ValueType SW, M, S;
	ValueType delta, r, SW_tmp;
};


#endif /* INCREMENTALWEIGHTEDMEANVARIANCECOMPUTER_H */
