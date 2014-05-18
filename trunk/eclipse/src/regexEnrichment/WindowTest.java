/**
 * 
 */
package regexEnrichment;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

/**
 * @author berald01
 *
 */
public class WindowTest {

	
	/**
	 * Test method for {@link regexEnrichment.Window#getChrom()}.
	 */
	@Test

	public void testContigWindows() {

		List<Window> windowCoords= new ArrayList<Window>();

		Window w= new Window();
		w.setStart(0);
		w.setEnd(15);
		w.setSequence("AACCTTGGNNaaaaa");
		System.out.println(w);
		windowCoords.add(w);

		w= new Window();
		w.setStart(10);
		w.setEnd(25);
		w.setSequence("AACCTTGGNNaaaaa");
		System.out.println(w);
		windowCoords.add(w);

		System.out.println(windowCoords.toString());

	}

}
