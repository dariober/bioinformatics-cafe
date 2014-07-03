package aligner;
import java.io.*;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import org.biojava3.alignment.*;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.core.sequence.*;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;

/**
 * MEMO: In each pair returned by Alignments.getAllPairsScores: 
 * The first in pair is called 
 *     QUERY by biojava. We call it QUERY or READ
 * The second in pair is called
 * 	   TARGET by biojava. We call it REFERENCE
 * @author berald01
 *
 */
public class Aligner {

	private SequencePair<DNASequence, NucleotideCompound> pair;
	private char strand;
	private int score;
	
	public SequencePair<DNASequence, NucleotideCompound> getPair() {
		return pair;
	}

	public void setPair(SequencePair<DNASequence, NucleotideCompound> pair) {
		this.pair = pair;
	}

	public char getStrand() {
		return strand;
	}

	public void setStrand(char strand) {
		this.strand = strand;
	}

	public int getScore() {
		return score;
	}

	public void setScore(int score) {
		this.score = score;
	}
	
	public void align(String query, String ref, String matrixFile) throws FileNotFoundException {
	
		AmbiguityDNACompoundSet ambiguityDNACompoundSet= new AmbiguityDNACompoundSet();
		
		DNASequence read= new DNASequence(query, ambiguityDNACompoundSet);
		DNASequence readRC= new DNASequence(read.getReverseComplement().getSequenceAsString(), ambiguityDNACompoundSet);

		DNASequence reference= new DNASequence(ref, ambiguityDNACompoundSet);
	
		SubstitutionMatrix<NucleotideCompound> matrix = new SimpleSubstitutionMatrix<NucleotideCompound>(ambiguityDNACompoundSet, 
				new File(matrixFile));
		
		// Forward
		List<DNASequence> lst= Arrays.asList(reference, read);
		int score= Alignments.getAllPairsScores(lst,
				PairwiseSequenceScorerType.GLOBAL, new SimpleGapPenalty(), matrix)[0];

		// Reverse complement
		List<DNASequence> lstRC= Arrays.asList(reference, readRC);
		int scoreRC= Alignments.getAllPairsScores(lstRC,
				PairwiseSequenceScorerType.GLOBAL, new SimpleGapPenalty(), matrix)[0];
		
		DNASequence bestRead;
		int bestScore;
		char strand;
		if (score > scoreRC){
			bestRead= read;
			bestScore= score;
			strand= '+';
		} else if (score < scoreRC){
			bestRead= readRC;
			bestScore= scoreRC;
			strand= '-';
		} else {
			Random rand= new Random(); 
			if (rand.nextBoolean()){
				bestRead= read;
				bestScore= score;
				strand= '+';
			} else {
				bestRead= readRC;
				bestScore= scoreRC;
				strand= '-';
			}
		}
		
		SimpleGapPenalty gap= new SimpleGapPenalty(); // Open default= -10; Extend default= -1 
		
		SequencePair<DNASequence, NucleotideCompound> pair= Alignments.getPairwiseAlignment(
				bestRead,   // query
				reference,  // target
				PairwiseSequenceAlignerType.GLOBAL, gap, matrix);
		
		this.pair= pair;
		this.strand= strand;
		this.score= bestScore;
	
	}

	
	
	public String getTmpFileFromResourceFile(String filename){
        
        //Open input file as stream
        BufferedReader br= null;
        String outfile= null;
        try {
                // Read script file
                br = new BufferedReader(
                                new InputStreamReader(getClass().getResourceAsStream(filename)));
       
                // Prepare tmp output
                File temp = File.createTempFile(filename, ".tmp");
                temp.deleteOnExit();
               
                // Write script to tmp
                FileWriter fw= new FileWriter(temp);
                String line;
                while ((line= br.readLine()) != null){
                        fw.write(line + "\n");
                }
                br.close();
                fw.close();
               
                // Return full path to tmp
                outfile= temp.getAbsolutePath();
               
        } catch (IOException e1) {
                e1.printStackTrace();
        }
        return(outfile);
	}

	public String toString(){

		StringBuilder sb= new StringBuilder();
		
		sb.append("Q: " + pair.getAlignedSequences().get(0).getSequenceAsString()); // query
		sb.append("\n");
		sb.append("R: " + pair.getAlignedSequences().get(1).getSequenceAsString()); // reference (or target)
		sb.append("\t");
		sb.append(score);
		sb.append("\t");
		sb.append(strand);
		sb.append("\t");
		sb.append(pair.getLength());
		
		return(sb.toString());
				
	}
		
}
