package bisReadBias;

import static org.junit.Assert.*;

import java.io.IOException;
import java.util.*;

import org.junit.*;

public class ScriptRunnerTest {

	@Test
	public void testSpeedList(){

		long t0= System.currentTimeMillis();
		int i= 0;
		int size= 10000;
		while(i < 1){
			List<Integer> list= new ArrayList<Integer>(size);
			int j= 0;
			while(j < size){
				list.add(j);
				j++;
			}
			i++;
		}
		long t1= System.currentTimeMillis();
		System.out.println("List: " + " " + (t1-t0));
	}
	
	@Test
	public void testSpeedArray(){

		long t0= System.currentTimeMillis();
		int i= 0;
		int size= 100000000;
		while(i < 1){
			int[] array= new int[size];
			int j= 0;
			while(j < size){
				array[j]= j;
				j++;
			}
			i++;
		}
		long t1= System.currentTimeMillis();
		System.out.println("Array: " + " " + (t1-t0));
	}
	
	
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
