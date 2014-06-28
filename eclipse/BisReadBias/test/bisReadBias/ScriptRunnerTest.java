package bisReadBias;

import static org.junit.Assert.*;

import java.io.IOException;

import org.junit.*;

public class ScriptRunnerTest {

	@Test
	public void givePathToResourceFile() throws IOException, InterruptedException{
		
		ScriptRunner sr= new ScriptRunner();
		String tmpFile= null;

		try {
			tmpFile= sr.getTmpFileFromResourceFile("Plotter.R");
		} catch(NullPointerException e){
			e.printStackTrace();
		}
		assertNotNull(tmpFile);

	}
	
	@Test
	public void canExecuteRscriptToPrintHelpToStderr() throws IOException, InterruptedException{
		
		ScriptRunner sr= new ScriptRunner();
		sr.runRscript("", "");
		assertEquals(0, sr.getReturnCode());
	
	}

	@Test
	public void canExecutePlotter() throws IOException, InterruptedException{
		
		ScriptRunner sr= new ScriptRunner();
		sr.runRscript(sr.getTmpFileFromResourceFile("Plotter.R"), "");
		assertEquals(100, sr.getReturnCode());
	
	}
	
	@Test
	public void readProfileFileAndPlot() throws IOException, InterruptedException{
		
		ScriptRunner sr= new ScriptRunner();
		sr.runRscript(sr.getTmpFileFromResourceFile("Plotter.R"), "test_data/profile_test.txt");
		assertEquals(0, sr.getReturnCode());
	
	}
	
}
