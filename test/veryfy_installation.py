from setuptools import setup, Command

class VerifyInstall(Command):
    description = 'Verify the installation of the package.'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        import subprocess
        subprocess.check_call(['python', 'tests/verify_installation.py'])

setup(
    name='MicrobeRX',
    version='1.0.0',
    packages=['microberx'],
    # Other arguments...
    cmdclass={
        'verify_install': VerifyInstall,
    },
)