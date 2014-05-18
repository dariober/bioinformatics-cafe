package markovChain;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

import markovChain.ArgParse;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

public class Main {

	public static void main(String[] args) {

		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		
		String fasta= opts.getString("fasta");
		String output= opts.getString("output");
		int nrandom= opts.getInt("nrandom");
		int norder= opts.getInt("norder");
		int length= opts.getInt("length");
		double pct_pseudocounts= opts.getDouble("pct_pseudocounts");
		/* End parsing command line arguments */

		/* Check and set command line arguments */
		// ...
		//
		/* ------------------------------------------------------------------ */

		File ff= new File(fasta);
		FastaSequenceFile fa= new FastaSequenceFile(ff, false);
		
		// Prepare output
		BufferedWriter wr= new BufferedWriter(Tools.printer(output));
		
		Alphabet alphabet= new Alphabet();
	
		while(true){
			// Iterate trough each sequence in the fasta file.
			ReferenceSequence faseq= fa.nextSequence();
			if (faseq == null){
				break;
			}
			String seq= new String(faseq.getBases()).toUpperCase();

			// Produce model for input sequence
			Set<Character> alphaset= alphabet.getAlphabetFromString(seq);
			Model model= new Model(norder, alphaset, seq, pct_pseudocounts);

			// Length of random sequence to generate
			int rndSeqLength= (length <= 0 ? seq.length() : length);
			
			MarkovSequence markovSequence= new MarkovSequence();
			int nseq= 1;
			while (nseq <= nrandom){
				// Generate chains
				String chain= markovSequence.generateMarkovSequence(model, rndSeqLength);
				try {
					wr.write(">" + faseq.getName() + ":" + nseq + "\n");
					wr.write(chain + "\n");
				} 
				catch(IOException e){
					e.printStackTrace();
				}
				nseq += 1;
			}
		} // End while loop through fasta seqs
		try{
			wr.close();
		}
		catch(IOException e){
			e.printStackTrace();
		}
	}
}
