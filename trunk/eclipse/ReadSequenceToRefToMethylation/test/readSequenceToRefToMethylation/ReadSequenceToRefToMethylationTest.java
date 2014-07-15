package readSequenceToRefToMethylation;

import static org.junit.Assert.*;

import org.junit.Test;

public class ReadSequenceToRefToMethylationTest {

	@Test
	public void nucleotideToMethylation() {
		AlignedRead aln= new AlignedRead();
		assertEquals('M', aln.nucleotideToMethylation('C', 'C', true, true));
		assertEquals('u', aln.nucleotideToMethylation('t', 'C', true, true));
		assertEquals('*', aln.nucleotideToMethylation('A', 'C', true, true));
		assertEquals('.', aln.nucleotideToMethylation('A', 'A', true, true));

		assertEquals('M', aln.nucleotideToMethylation('C', 'C', false, false));
		assertEquals('u', aln.nucleotideToMethylation('t', 'C', false, false));
		assertEquals('*', aln.nucleotideToMethylation('A', 'C', false, false));
		assertEquals('.', aln.nucleotideToMethylation('A', 'A', false, false));
	
		assertEquals('M', aln.nucleotideToMethylation('G', 'G', false, true));
		assertEquals('u', aln.nucleotideToMethylation('A', 'G', false, true));
		assertEquals('*', aln.nucleotideToMethylation('C', 'G', true, false));
		assertEquals('*', aln.nucleotideToMethylation('C', 'G', false, true));

	}

}
