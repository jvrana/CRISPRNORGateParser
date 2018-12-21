from distutils.core import setup

setup(name='NORGateParser',
      packages=['norgateparser'],
      version='0.0.1',
        package_data={
            'norgateparser': {
                'data/*'
            }
        }
     )
