import os
import subprocess
import webbrowser
import threading
import time
import socket
import sys
from pathlib import Path
import shutil

def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        s.listen(1)
        port = s.getsockname()[1]
    return port

def check_node_npm():
    """Check if Node.js and npm are available."""
    try:
        subprocess.run(['node', '--version'], check=True, capture_output=True)
        subprocess.run(['npm', '--version'], check=True, capture_output=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def install_dependencies(gui_dir):
    """Install npm dependencies if node_modules doesn't exist."""
    node_modules = gui_dir / "node_modules"
    if not node_modules.exists():
        print("Installing GUI dependencies...")
        result = subprocess.run(['npm', 'install'], cwd=gui_dir, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Failed to install dependencies: {result.stderr}")
        print("Dependencies installed successfully!")

def build_gui(gui_dir):
    """Build the Next.js application."""
    print("Building GUI (this may take a few minutes on first run)...")
    result = subprocess.run(['npm', 'run', 'build'], cwd=gui_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print("Build failed. Error output:")
        print(result.stderr)
        raise RuntimeError(f"Failed to build GUI: {result.stderr}")
    print("GUI built successfully!")

def launch_gui(port=None, open_browser=True, dev_mode=False):
    """Launch the SIMBA Chem GUI."""
    if not check_node_npm():
        raise RuntimeError(
            "Node.js and npm are required to run the GUI. "
            "Please install Node.js from https://nodejs.org/"
        )
    
    if port is None:
        port = find_free_port()
    
    gui_dir = Path(__file__).parent / "gui"
    if not gui_dir.exists():
        raise FileNotFoundError(
            "GUI files not found. Make sure the package was installed correctly."
        )
    
    # Install dependencies if needed
    install_dependencies(gui_dir)
    
    # Check if we need to build (only for production mode)
    if not dev_mode:
        next_dir = gui_dir / ".next"
        if not next_dir.exists():
            build_gui(gui_dir)
    
    # Set up environment
    env = os.environ.copy()
    env['PORT'] = str(port)
    env['NODE_ENV'] = 'development' if dev_mode else 'production'
    
    # Choose command based on mode
    if dev_mode:
        cmd = ['npm', 'run', 'dev']
        startup_time = 5  # Dev server takes longer to start
        print(f"Starting SIMBA Chem GUI in development mode on port {port}...")
    else:
        cmd = ['npm', 'start']
        startup_time = 3
        print(f"Starting SIMBA Chem GUI on port {port}...")
    
    # Start Next.js server
    process = subprocess.Popen(
        cmd,
        cwd=str(gui_dir),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        universal_newlines=True
    )
    
    # Wait for server to start and optionally open browser
    if open_browser:
        def open_browser_delayed():
            time.sleep(startup_time)
            if process.poll() is None:  # Check if process is still running
                webbrowser.open(f'http://localhost:{port}')
        
        browser_thread = threading.Thread(target=open_browser_delayed)
        browser_thread.start()
    
    print(f"GUI will be available at http://localhost:{port}")
    print("Waiting for server to start...")
    
    # Show initial output to indicate progress
    server_ready = False
    try:
        for line in process.stdout:
            print(f"GUI: {line.strip()}")
            
            # Check for server ready indicators
            if any(indicator in line.lower() for indicator in ['ready', 'compiled', 'started server']):
                if not server_ready:
                    print(f"\n✓ GUI server is ready at http://localhost:{port}")
                    print("Press Ctrl+C to stop the server\n")
                    server_ready = True
        
        process.wait()
        
    except KeyboardInterrupt:
        print("\nShutting down GUI server...")
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            print("Force killing server...")
            process.kill()
        print("GUI server stopped.")
    except Exception as e:
        print(f"Error running GUI: {e}")
        process.terminate()
    finally:
        if process.poll() is None:
            process.terminate()

def launch_gui_dev():
    """Launch GUI in development mode with hot reloading."""
    launch_gui(dev_mode=True)

if __name__ == "__main__":
    launch_gui()