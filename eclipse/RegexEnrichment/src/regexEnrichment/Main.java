package regexEnrichment;

import net.sf.picard.reference.*;
import net.sourceforge.argparse4j.inf.Namespace;

import com.google.common.collect.*;

import java.io.*;
import java.util.*;
import java.util.regex.*;

public class Main {
	public static void main(String[] args){
			
		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		
		String fasta= opts.getString("fasta");
		String output= opts.getString("output");
		int window_size= opts.getInt("window_size");
		String pattern= opts.getString("pattern");
		int orderMarkovChain= opts.getInt("orderMarkovChain");
		double pctPseudocounts= opts.getDouble("pctPseudocounts");
		boolean casesensitive= opts.getBoolean("casesensitive");
		int precision= opts.getInt("precision");
		int nrand= opts.getInt("nrand");
		boolean revcomp= !opts.getBoolean("norevcomp");
		String seqcase= opts.getString("seqcase");
		int step= ( (opts.getInt("step") <= 0) ? window_size : opts.getInt("step") );
		String wLookUpTable= opts.getString("wLookUpTable");
		String rLookUpTable= opts.getString("rLookUpTable");
		boolean verbose= opts.getBoolean("verbose");
		
		if (!casesensitive){
			seqcase= "U";
		}
		/* End parsing arguments */

		// Prepare output
		BufferedWriter wr= new BufferedWriter(Tools.printer(output));
		// ---
		
		File ff= new File(fasta);
		FastaSequenceFile fa= new FastaSequenceFile(ff, false); 
		
		// MEMO: Flag 2 means CASE_INSENSITIVE.
		int isCaseSens= (casesensitive ? 0 : 2);
		Pattern cpattern= Pattern.compile(pattern, isCaseSens); 
				
		ProbsLookUpTable ptable= new ProbsLookUpTable();
		if (rLookUpTable != null){
			ptable= ptable.readMap(rLookUpTable);
			if (ptable.getNrand() < nrand){
				System.err.println("Number of required randomizations is greater than in look-up file");
				System.exit(1);
			}
			if (!ptable.getPattern().equals(pattern)){
				System.err.println("Pattern does not equal the pattern in look-up file");
				System.err.println("Required: " + pattern);
				System.err.println("Look up: " + ptable.getPattern());
				System.exit(1);
			}
			if (ptable.getWindow_size() != window_size){
				System.err.println("Window size does not equal the number in look-up file");
				System.exit(1);
			}
			if (ptable.getPrecision() != precision){
				System.out.println(ptable);
				System.err.println("Precision different than in look-up file");
				System.err.println("Look up: " + ptable.getPrecision());
				System.err.println("Required: " + precision);
				System.exit(1);
			}
		}
		else {
			ptable.setNrand(nrand);
			ptable.setPattern(pattern);
			ptable.setWindow_size(window_size);
			ptable.setPrecision(precision);
		}
		ContigWindows contigWindows= new ContigWindows(); 
				
		int i= 0;
		long t1;
		long t0= System.currentTimeMillis();
		long tz= System.currentTimeMillis();
		while(true){
			// Iterate trough each sequence in the fasta file.
			ReferenceSequence seq= fa.nextSequence();
			if (seq == null){
				break;
			}
			// Set sequence case
			String faseq= null;
			if (seqcase.equals("U")){
				faseq= new String(seq.getBases()).toUpperCase();
			}
			else if (seqcase.equals("L")){
				faseq= new String(seq.getBases()).toLowerCase();
			}
			else if (seqcase.equals("k")){
				faseq= new String(seq.getBases());
			} else {
				System.err.println("Invalid case");
				System.exit(1);
			}
			
			// Divide entire sequence in windows of required size and overlap:
			contigWindows= new ContigWindows(seq.getName(), faseq);
			contigWindows.setWindowCoords(window_size, step);
			for(Window w : contigWindows.getWindowCoords()){
				// Iterate through each window
				i++;
				w.setNucleotideFreq(precision);
				List<SeqMatcher> matches= Tools.regexMatcher(cpattern, w.getSequence(), revcomp);
				/* Decide if randomization has to run */						

				boolean moreHitsThanEver= ptable.obsHitsGreaterThanMax(matches.size(), w.getNucleotideFreq());
				
				// Note: If orderMarkoVModel is >0, pvalues are always recomputed even if the 
				// nucleotide composition is the same. The Markov model depends on the composition
				// and the order.
				if (orderMarkovChain > 0 || !ptable.probsLookUpTable.containsKey(w.getNucleotideFreq()) || moreHitsThanEver){
					TreeMultiset<Integer> nullDistr= TreeMultiset.create();
					if (ptable.probsLookUpTable.containsKey(w.getNucleotideFreq()) ){
						// Get collection of hits from previous loop.
						// This can happen if the current nucleotide composition has been encountered before
						// but not enough randomizations where done for the current number of actual hits.
						nullDistr= ptable.probsLookUpTable.get(w.getNucleotideFreq());
					}
					NullDistribution nullDistribution= new NullDistribution(
							w.getSequence(), orderMarkovChain, pctPseudocounts);
					nullDistr= nullDistribution.generateNullDistribution(
							w.getSequence(), 
							cpattern, 
							nrand, 
							precision, 
							matches.size(),
							nullDistr,
							revcomp
					);
					ptable.probsLookUpTable.put(w.getNucleotideFreq(), nullDistr);
					ptable.maxHitsLookUpTable.put(w.getNucleotideFreq(), matches.size());
				}
				TreeMultiset<Integer> nullHits= ptable.probsLookUpTable.get(w.getNucleotideFreq());
				float pvalue= Tools.mapObservedHitsToNull(nullHits, matches.size());
				String xpvalue= Tools.formatPvalue(pvalue, nullHits.size());
				for (SeqMatcher m : matches ){
					// Print results
					String bed= new BedLine(w, m, xpvalue).toString();
					try {
						t1= System.currentTimeMillis();
						wr.write(bed + "\t" + nullHits + "\n");
						wr.flush();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}		
				t1= System.currentTimeMillis();
				if (verbose && (t1 - tz) > 60000){
					Tools.pregressReport(t0, t1, i, w, ptable);
					tz= t1;
				}

			} // End of this sequence loop
		} // End Fasta sequences loop
		if (wLookUpTable != null){
			ptable.writeMap(wLookUpTable);
		} 
		System.exit(0);
	} // end main methid
} // end Main class
