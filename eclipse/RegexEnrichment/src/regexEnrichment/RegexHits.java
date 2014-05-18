package regexEnrichment;

import java.util.*;
import java.util.regex.*;

import com.google.common.collect.*;

public class RegexHits {
	private int observedRegexHits;
	private HashMultiset<Integer> nullDistr; 
	private float pvalue;

	public int getObservedRegexHits() {
		return observedRegexHits;
	}
	public void setObservedRegexHits(int observedRegexHits) {
		this.observedRegexHits = observedRegexHits;
	}
	public HashMultiset<Integer> getNullDistr() {
		return nullDistr;
	}
	public void setNullDistr(HashMultiset<Integer> nullDistr) {
		this.nullDistr = nullDistr;
	}
	public float getPvalue() {
		return pvalue;
	}
	public void setPvalue(float pvalue) {
		this.pvalue = pvalue;
	}
	
	public float pvalOfObsHits(Pattern pattern, Window window, boolean revcomp, int precision, 
			Map<Map<String, Float>, TreeMultiset<Integer>> probLookUpMap){
//		List<SeqMatcher> matches= Tools.regexMatcher(pattern, window.getSequence(), revcomp);
		observedRegexHits= Tools.regexCounter(pattern, window.getSequence(), revcomp);
		float pvalue= Tools.mapObservedHitsToNull(probLookUpMap.get(window.getNucleotideFreq()), observedRegexHits);
		return(pvalue);
	}
}
