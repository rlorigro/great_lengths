import pyfastx
import os


def main():
    project_directory = os.path.abspath(os.path.dirname(__file__))
    test_directory = os.path.join(project_directory, "data/test")
    test_filename = "simple.fastq"
    test_path = os.path.join(test_directory, test_filename)

    print("project_directory:", project_directory)
    print("test_directory:", test_directory)
    print("test_path:", test_path)

    fastq = pyfastx.Fastq(test_path)

    for s in fastq:
        print(s.len)


if __name__ == "__main__":
    main()
