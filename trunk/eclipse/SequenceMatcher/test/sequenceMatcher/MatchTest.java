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
		
		matches.setHD(0);
		assertTrue(matches.getHD() == -1);
		
		matches.setLD(1);
		assertTrue(matches.getLD() == 1);
		
		matches.setJWD();
		assertTrue(matches.getJWD() == 0.92);
		
		double d= matches.getFilterDistance("HD", -1);
		assertTrue(d == 1);

		d= matches.getFilterDistance("LD", 1);
		assertTrue(d == 1);
		
		d= matches.getFilterDistance("LD", -1);
		assertTrue(d == 1);
		
		d= matches.getFilterDistance("LD", 0);
		assertTrue(d == -1);
		
		matches= new Match(seqA, seqA);
		d= matches.getFilterDistance("LD", 0);
		assertTrue(d == 0);

		d= matches.getFilterDistance("LD", 1);
		assertTrue(d == 0);

		d= matches.getFilterDistance("LD", -1);
		assertTrue(d == 0);
		
	}
}
