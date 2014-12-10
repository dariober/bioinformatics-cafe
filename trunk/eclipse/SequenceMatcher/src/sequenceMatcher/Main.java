package sequenceMatcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;

import org.biojava3.core.sequence.DNASequence;

import net.sourceforge.argparse4j.inf.Namespace;

public class Main {

	public static void main(String[] args) throws IOException {
		
		/* Start parsing arguments */
		Namespace opts= ArgParse.argParse(args);
		ArgParse.validateArgs(opts);
		String a= opts.getString("a");
		String b= opts.getString("b");
		String method= opts.getString("method");
		int nm= opts.getInt("nm");
		boolean norc= opts.getBoolean("norc");
		
		/* Adjust args */
		
		// How many loops thorough each sequence in B? 
		// 1 if rev comp is not required. 2 otherwise (+ and -)
		int nLoops= (!norc) ? 2 : 1;
		
//		int NMAX= 0; // Max edit distance to return a match
//		String faA= "/Users/berald01/Tritume/setA.fa"; // All read in memory
//		String faB= "/Users/berald01/Tritume/setAB.fa"; // Read one sequence at a time
//		String method= "Hamming";
//		boolean revcomp= true;
		// ---------------------------------------------------------------------
	
		ArrayList<String[]> fastaFileA= SequenceReader.readFastaToList(a);
				
		System.err.println("Sequences in A: " + fastaFileA.size());
		
		// ---------------------------------------------------------------------
		// Prepare reading file b
		BufferedReader brB= Opener.openBr(b);
				
		long t0= System.currentTimeMillis();
		
		System.out.println(new Match().getHeader());
		
		int nseq= 0;
		int nmatch= 0;
		String[] seqB;
		while((seqB= SequenceReader.getNextSequence(brB)) != null){

			String seqBrc= new DNASequence(seqB[1]).getReverseComplement().getSequenceAsString();
			
			for(int i=0; i < fastaFileA.size(); i++){
				int nLoopCnt= 0;
				while(nLoopCnt < nLoops) {
					String[] seqA= fastaFileA.get(i);
					Match m= new Match(seqA, seqB);
					if(nLoopCnt == 0){
						m.setStrand("+");
					} else if (nLoopCnt == 1){
						m.setStrand("-");
						m.setSeqB(seqBrc);
					} else {
						System.exit(1);
					}
					m.getDistance(method, nm);
					if(nm < 0 || (m.getNM() >= 0 && m.getNM() <= nm)){
						m.align();
						nmatch++;
						System.out.println(m);
					}
					nLoopCnt++;
				}
			} // end for loop file A
			nseq++;
			if(nseq % 1000 == 0){
				System.err.println("Processed " + nseq + " sequences from B");
			}
		}
		
		long t1= System.currentTimeMillis();
		System.err.println("Completed in " + (t1-t0)/1000.0 + "s");
		System.err.println("N matches: " + nmatch);		

	}

}


//String[] seqA= fastaFileA.get(i);
//Match m= new Match(seqA, seqB);
//m.getDistance(method, nm);
//if(nm < 0 || (m.getNM() >= 0 && m.getNM() <= nm)){
//	nmatch++;
//	m.align();
//	System.out.println(m);
//}
//if(!norc){
//	// Replace sequence A with its revcomp.
//	m.setSeqA(seqA[2]);
//	m.getDistance(method, nm);
//	m.setStrand("-");
//	if(nm < 0 || (m.getNM() >= 0 && m.getNM() <= nm)){
//		nmatch++;
//		System.out.println(m);
//	}
//}
