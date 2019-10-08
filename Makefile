install:
	# Install the python script to the bin
	cp precip_rescaling.py /usr/local/bin/precip_rescaling
	chmod +x /usr/local/bin/precip_rescaling

uninstall:
	# Remove the python script from the bin
	rm /usr/local/bin/precip_rescaling

develop:
	# Link the python script so edits are reflected in realtime
	# make uninstall
	chmod +x precip_rescaling.py
	ln precip_rescaling.py /usr/local/bin/precip_rescaling
