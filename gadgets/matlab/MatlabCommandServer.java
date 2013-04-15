import java.lang.Thread;
import java.net.ServerSocket;
import java.net.Socket;
import java.io.*;
import com.mathworks.jmi.*;

class MatlabCommandServer extends Thread {

    ServerSocket socket = null;
    Matlab matlab = null;
    boolean stop_signal_received = false;
    boolean socket_is_open = false;

    public MatlabCommandServer() {
	try {
	    matlab = new Matlab();
	} catch (Exception e) {
	    System.err.println("Failed to create new Matlab");
            System.err.println(e.getMessage());
	}
    }

    private boolean openSocket() {
        if (socket_is_open) {
            return false;
        }

        try {
            socket = new ServerSocket(3000);
            socket.setSoTimeout(1000); //1000ms time out. We will check if shutdown has occurred every 1000ms.
        } catch (Exception e) {
            // Socket creation has failed, we should do something
            System.err.println("Socket failed to open");
            System.err.println(e.getMessage());
            return false;
        }
        socket_is_open = true;
        return true;
    }

 
    private boolean closeSocket() {
        if (!socket_is_open) {
            return false;
        }

        try {
            socket.close();
        } catch (Exception e) {
            // Socket close has failed, we should do something
            System.err.println("Socket failed to close");	    
            System.err.println(e.getMessage());
            return false;
        }
        socket_is_open = false;
        return true;
    }

    private boolean receiveCommand() {
	try {
	    Socket sock = socket.accept();
	    BufferedReader in = new BufferedReader(new InputStreamReader(sock.getInputStream()));

	    //System.out.println("Waiting for command");
	    while (!in.ready()) ;

	    String command = in.readLine();
        
	    //System.out.println(command);
	    matlab.evalConsoleOutput(command);

	    in.close();
	    sock.close();


	} catch (java.io.InterruptedIOException e) {
             // This means that we have waited for connection but so far nothing. 
             // We should check if the thread has been notified to stop,
             // if so, stop the loop and otherwise continue.
	    if (stop_signal_received) {
		return false;
	    }
	} catch (Exception e) {
	    System.err.println("Something unexpected has happened!!");
	    System.err.println(e.getMessage());
	    return false;
	}
	return true;

    }

    public void notifyEnd() {
	stop_signal_received = true;
    }

    public void run() {
	if (!openSocket()) {
            return;
        }

	while (true) {
	    if (!receiveCommand()) break;
	}
	closeSocket();
	stop_signal_received = false;
    }

    protected void finalize() throws Throwable {
	System.out.println("MatlabMessageServer finalize() called");
	stop_signal_received  = true;
	closeSocket();
	super.finalize();
    }

}
