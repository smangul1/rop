from subprocess import call, Popen, PIPE
checksum_one = "8cafeff914e8448e5659d817c294abc3"
checksum_test = Popen(["md5sum", "install.py"], stdin=PIPE)
print checksum_test