package sequenceMatcher;

import static org.junit.Assert.*;

import org.junit.Test;

public class MatchTest {

	@Test
	public void canGetDistance() {
		//Create a set of sequences and a single sequence. 
		String[] seqA= {"seq1-A", "ACTGA"};
		String[] seqB= {"seq1-B", "ACTGN"};

		Match matches= new Match(seqA, seqB);
		
		double jwd= matches.getJaroWinklerDistance();
		System.out.println(jwd);
		
		matches.getDistance("Leven", -1);
		assertTrue(matches.getNM() == 1);
		
		matches.getDistance("Leven", 0);
		assertTrue(matches.getNM() == -1);

		matches.getDistance("Hamming", -1);
		assertTrue(matches.getNM() == 1);
		
		matches.getDistance("Hamming", 0);
		assertTrue(matches.getNM() == -1);
		
	}
}
