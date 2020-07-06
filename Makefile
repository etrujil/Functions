install:
	# Install the python script to the bin
	cp precip_rescaling.py /usr/local/bin/precip_rescaling
	chmod +x /usr/local/bin/precip_rescaling

	cp nc_replace.py /usr/local/bin/nc_replace
	chmod +x /usr/local/bin/nc_replace

	cp image_calc.py /usr/local/bin/image_calc
	chmod +x /usr/local/bin/image_calc


uninstall:
	# Remove the python script from the bin
	rm /usr/local/bin/precip_rescaling
	rm /usr/local/bin/nc_replace
	rm /usr/local/bin/image_calc


develop:
	# Link the python script so edits are reflected in realtime
	chmod +x precip_rescaling.py
	ln precip_rescaling.py /usr/local/bin/precip_rescaling

	chmod +x nc_replace.py
	ln nc_replace.py /usr/local/bin/nc_replace

	chmod +x image_calc.py
	ln image_calc.py /usr/local/bin/image_calc