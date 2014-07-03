package aligner;

import static org.junit.Assert.*;

import java.io.FileNotFoundException;

import org.biojava3.core.sequence.DNASequence;
import org.junit.Test;

public class AlignerTest {
	
	@Test
	public void alignToReverseNoGaps() throws FileNotFoundException{
		
		String query= "ATCGAGAATCCCGGTGCCGATACCNACTCTTGNAGAA";
		String ref=   "TTCTNCAAGAGTNGGTATCGGCACCGGGATTCTCGAT";

		Aligner aln= new Aligner();
		String matrixFile= aln.getTmpFileFromResourceFile("nuc_N1.txt");
		
		aln.align(query, ref, matrixFile);
		
		assertEquals('-', aln.getStrand());
		
		//First in pair is query. revcomp'd to match reference
		assertEquals(new DNASequence(query).getReverseComplement().getSequenceAsString(), aln.getPair().getAlignedSequences().get(0).getSequenceAsString());

		// Second in pair is reference.
		assertEquals(ref, aln.getPair().getAlignedSequences().get(1).getSequenceAsString());
		
	}
	
	@Test
	public void alignUnequalLength() throws FileNotFoundException{
		
		String query= "GGGGGGGGGGATCGAGAATCCCGGTGCCGATACCNACTCTTGNAGAATTTTTT";
		String ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCNACTCTTGNAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";

		Aligner aln= new Aligner();
		String matrixFile= aln.getTmpFileFromResourceFile("nuc_N1.txt");
		
		aln.align(query, ref, matrixFile);
		
		assertEquals('+', aln.getStrand());
		
		System.out.println(aln.getPair().getAlignedSequences().get(0).getSequenceAsString());
		System.out.println(aln.getPair().getAlignedSequences().get(1).getSequenceAsString());
		
	}

	@Test
	public void targetIsReference() throws FileNotFoundException{

		String query= "AAAhCCC";
		String ref=   "AAAnCCC";
		Aligner aln= new Aligner();
		String matrixFile= aln.getTmpFileFromResourceFile("nuc_N1.txt");
		
		aln.align(query, ref, matrixFile);
		
		String referenceBase= aln.getPair().getCompoundInTargetAt(4).toString();
		assertEquals("n", referenceBase);
		
		String queryBase= aln.getPair().getCompoundInQueryAt(4).toString();
		assertEquals("h", queryBase);
				
	}

	
	@Test
	public void alignFullLength() throws FileNotFoundException{

		String query= "GGGAAGCAAAATCGAGAATCCCGGTGCCGATACCCACTCTTGCAGAATCTATCACAACGTGTCAGATATATACATCCGAT";
		String ref=   "NNNNNNNNNNATCGAGAATCCCGGTGCCGATACCYACTCTTGYAGAANNNNNNNNNNCGTGTCAGATATATACATCCGAT";

		Aligner aln= new Aligner();
		String matrixFile= aln.getTmpFileFromResourceFile("nuc_N0.txt");
		
		aln.align(query, ref, matrixFile);
		
		System.out.println(aln);

	}
}
