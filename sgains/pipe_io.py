import queue

from subprocess import Popen, PIPE
from threading import Thread


class InputGenerateThread(Thread):

    def __init__(self, control_queue, input_fd, input_function_generator):
        super(InputGenerateThread, self).__init__()

        self.control_queue = control_queue
        self.input_fd = input_fd
        self.input_function_generator = input_function_generator

    def run(self):
        for line in self.input_function_generator():
            self.input_fd.write(bytes(line, 'utf8'))
            self.input_fd.write(bytes("\n", 'utf8'))
            self.input_fd.flush()
        self.input_fd.close()
        self.control_queue.put("in_done")


class OutputProcessThread(Thread):

    def __init__(self, control_queue, output_fd, output_process_function):
        super(OutputProcessThread, self).__init__()

        self.control_queue = control_queue
        self.output_fd = output_fd
        self.output_process_function = output_process_function

    def run(self):
        while True:
            line = self.output_fd.readline()
            if not line:
                print("output fd closed...")
                self.output_fd.close()
                self.control_queue.put("out_done")
                return

            self.output_process_function(line)


class PipeIO(object):

    def __init__(
            self, command,
            input_generate_function, output_process_function):
        self.command = command
        self.input_generate_function = input_generate_function
        self.output_process_function = output_process_function
        self.process = None
        self.control_queue = queue.Queue()

    def run(self):
        with Popen(self.command, stdout=PIPE, stdin=PIPE) as proc:

            self.input_thread = InputGenerateThread(
                self.control_queue,
                proc.stdin,
                self.input_generate_function)
            self.output_thread = OutputProcessThread(
                self.control_queue,
                proc.stdout,
                self.output_process_function)

            self.input_thread.start()
            self.output_thread.start()

            while True:
                msg = None
                try:
                    msg = self.control_queue.get()
                except queue.Empty:
                    print("timeout - queue empty")
                    msg = None
                if msg == 'out_done':
                    print("output done")
                    break
                if msg == 'in_done':
                    print('input done')
            self.input_thread.join()
            self.output_thread.join()


if __name__ == "__main__":
    def test_input_function():
        return ["line1", "line2"]

    def test_output_function(line):
        assert line is not None
        print("processing line:", line)

    test_pipe = PipeIO(
        ['cat'],
        test_input_function,
        test_output_function
    )
    test_pipe.run()
